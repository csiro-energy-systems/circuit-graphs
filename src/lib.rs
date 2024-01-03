pub mod circuit_graph {
    use std::collections::HashMap;

    use faer::prelude::SpSolver;
    use faer::solvers::Svd;
    use faer::{Col, Mat};
    use faer_entity::ComplexField;
    use petgraph::algo;
    use petgraph::prelude::*;

    #[derive(PartialEq, Eq)]
    pub enum VertexType {
        Source,
        Sink,
        Internal,
    }

    pub struct VertexMetadata<T> {
        pub voltage: Option<T>,
        pub tag: u32,
        pub vertex_type: VertexType,
        power: Option<T>,
    }

    impl<T> VertexMetadata<T> {
        /// Construct a new [`VertexMetadata`]. This constructor guarantees that a
        /// source or sink node will have a voltage associated, and that internal nodes
        /// will have [`None`] instead.
        pub fn new(voltage: Option<T>, tag: u32, vertex_type: VertexType) -> Self {
            match vertex_type {
                VertexType::Internal => Self {
                    voltage: None,
                    tag,
                    vertex_type,
                    power: None,
                },
                _ => {
                    assert!(!voltage.is_none());
                    Self {
                        voltage,
                        tag,
                        vertex_type,
                        power: None,
                    }
                }
            }
        }
    }

    /// A struct to hold the informatino for an edge in the graph.
    ///
    /// An edge represents a component in a circuit, which may be purely
    /// resistive, capacitative, or inductive, or some combination of the three.
    /// In any case, the admittance must be provided as a single value, whether
    /// that is a float (or some other real value) or a complex type. Complex
    /// values should be provided as `G + jB` where
    /// - `G` is conductance (real-valued),
    /// - `B` is susceptance (real-valued), and
    /// - `j` is the imaginary unit.
    ///
    /// In particular, the `c64` or `c32` types from `faer` are recommended for
    /// best performance.
    pub struct EdgeMetadata<T> {
        pub tail: u32,
        pub head: u32,
        pub admittance: T,
        current_id: Option<usize>,
        pub current: Option<T>,
        power: Option<T>,
    }

    impl<T> EdgeMetadata<T> {
        /// Constructs a new [`EdgeMetadata`] with the provided tail and head indices
        /// and admittance.
        ///
        /// See the [`EdgeMetadata`] documentation for details on admittance values.
        pub fn new(tail: u32, head: u32, admittance: T) -> Self {
            Self {
                tail,
                head,
                admittance,
                current_id: None,
                current: None,
                power: None,
            }
        }
    }

    /// A struct representing a passive circuit. It relies on a flow graph where
    ///  - vertices are annotated with a [`VertexType`],
    ///  - vertices are annotated with a voltage, and
    ///  - edges are annotated with a weight equal to the admittance (inverse of
    /// impedance) of the component they represent.
    ///
    /// While the use of 'admittance' implies a complex value, using a real type
    /// will work also for DC circuits or, in theory, purely resistive circuits.
    /// Complex admittance values should be in `G + jB` format; see the
    /// [`EdgeMetadata`] documentation for details.
    pub struct Circuit<T> {
        pub graph: DiGraph<VertexMetadata<T>, EdgeMetadata<T>>,
    }

    impl<T> Circuit<T> {
        pub fn new(vertices: Vec<VertexMetadata<T>>, edges: Vec<EdgeMetadata<T>>) -> Self {
            let mut graph: DiGraph<VertexMetadata<T>, EdgeMetadata<T>> = DiGraph::new();

            let mut vertex_indices: HashMap<u32, NodeIndex> = HashMap::new();

            for vertex in vertices {
                let tag = vertex.tag;
                let node_index = graph.add_node(vertex);
                vertex_indices.insert(tag, node_index);
            }

            for edge in edges {
                graph.add_edge(
                    *vertex_indices.get(&edge.tail).unwrap(),
                    *vertex_indices.get(&edge.head).unwrap(),
                    edge,
                );
            }

            Self { graph }
        }

        /// Find and label sets of edges which are in series. Returns the number of
        /// distinct currents which need to be found in the circuit.
        pub fn determine_unknown_currents(&mut self) -> usize {
            // This method uses the fact that every edge of the graph corresponds to a
            // component whose current we want to find, but that components in series (which
            // can be uniquely associated with vertices having deg_in = deg_out = 1)
            // need to be uncounted since any group all share the same current. Such nodes
            // are found here:
            let series_nodes: Vec<NodeIndex> = self
                .graph
                .node_indices()
                .filter(|index| {
                    self.graph.edges_directed(*index, Incoming).count() == 1
                        && self.graph.edges_directed(*index, Outgoing).count() == 1
                })
                .collect();

            let mut next_current_index: usize = 0;

            // Match up series nodes to all have the same current indices.
            for node_index in &series_nodes {
                // Note we know these edges exist since we just grabbed them from the graph, so
                // all the unwrap calls in this for block are safe.
                let incoming_edge_index = self
                    .graph
                    .edges_directed(*node_index, Incoming)
                    .last()
                    .unwrap()
                    .id();
                let outgoing_edge_index = self
                    .graph
                    .edges_directed(*node_index, Outgoing)
                    .last()
                    .unwrap()
                    .id();

                if let Some(id) = self
                    .graph
                    .edge_weight(incoming_edge_index)
                    .unwrap()
                    .current_id
                {
                    let outgoing_weight = self.graph.edge_weight_mut(outgoing_edge_index).unwrap();
                    outgoing_weight.current_id = Some(id);
                } else if let Some(id) = self
                    .graph
                    .edge_weight(outgoing_edge_index)
                    .unwrap()
                    .current_id
                {
                    let incoming_weight = self.graph.edge_weight_mut(incoming_edge_index).unwrap();
                    incoming_weight.current_id = Some(id);
                } else {
                    let incoming_weight = self.graph.edge_weight_mut(incoming_edge_index).unwrap();
                    incoming_weight.current_id = Some(next_current_index);

                    let outgoing_weight = self.graph.edge_weight_mut(outgoing_edge_index).unwrap();
                    outgoing_weight.current_id = Some(next_current_index);

                    next_current_index += 1;
                }
            }

            // Now that all the series pairs have been added, label every other edge in
            // increasing order.
            for edge_weight in self.graph.edge_weights_mut() {
                if edge_weight.current_id.is_none() {
                    edge_weight.current_id = Some(next_current_index);
                    next_current_index += 1;
                }
            }

            self.graph.edge_count() - series_nodes.len()
        }

        /// Count the number of unknown voltages that need to be found when solving.
        pub fn count_unknown_voltages(&self) -> usize {
            // Here we simply count the number of vertices with vertex_type Internal; sinks
            // and sources have their voltage predetermined.
            self.graph
                .node_weights()
                .filter(|v| v.vertex_type == VertexType::Internal)
                .count()
        }

        /// Find the sequences of edges which form source->sink paths within the circuit.
        pub fn find_paths(&self) -> Vec<(NodeIndex, Vec<EdgeIndex>)> {
            let source_indices: Vec<NodeIndex> = self.graph.externals(Incoming).collect();
            let sink_indices: Vec<NodeIndex> = self.graph.externals(Outgoing).collect();

            let mut paths: Vec<(NodeIndex, Vec<EdgeIndex>)> = Vec::new();

            for source_index in &source_indices {
                for sink_index in &sink_indices {
                    paths.append(&mut self.find_paths_between(&source_index, &sink_index));
                }
            }

            paths
        }

        /// Finds, as sequences of edges, all the paths between a given source and sink node
        /// using a modified depth-first search.
        /// Also returns the source node as a convenience.
        fn find_paths_between(
            &self,
            source_index: &NodeIndex,
            sink_index: &NodeIndex,
        ) -> Vec<(NodeIndex, Vec<EdgeIndex>)> {
            let mut paths: Vec<(NodeIndex, Vec<EdgeIndex>)> = Vec::new();
            let mut visited_edges: Vec<EdgeIndex> = Vec::new();
            let mut visited_nodes: Vec<NodeIndex> = Vec::new();

            let mut stack: Vec<_> = vec![self.graph.edges_directed(*source_index, Outgoing)];

            if source_index == sink_index {
                return paths;
            }

            while !stack.is_empty() {
                // Will never error since we know stack isn't empty
                let edge_iter = stack.last_mut().unwrap();
                match edge_iter.next() {
                    Some(next_edge) => {
                        if next_edge.target() == *sink_index {
                            let mut new_path = visited_edges.clone();
                            new_path.push(next_edge.id());
                            paths.push((source_index.clone(), new_path));
                        } else if !visited_nodes.contains(&next_edge.target()) {
                            visited_edges.push(next_edge.id());
                            visited_nodes.push(next_edge.source());
                            stack.push(self.graph.edges_directed(next_edge.target(), Outgoing));
                        }
                    }
                    None => {
                        stack.pop();
                        visited_edges.pop();
                        visited_nodes.pop();
                    }
                }
            }

            paths
        }
    }

    impl<
            T: ComplexField
                + std::ops::Add<Output = T>
                + std::ops::Sub<Output = T>
                + std::ops::Mul<Output = T>,
        > Circuit<T>
    {
        /// Solve current values
        pub fn solve_currents(&mut self) {
            let num_unknowns = self.determine_unknown_currents();

            let paths = self.find_paths();
            let num_paths = paths.len();

            // Now find internal vertices which bus multiple branches
            let bus_indices: Vec<NodeIndex> = self
                .graph
                .node_indices()
                .filter(|index| {
                    self.graph.node_weight(*index).unwrap().vertex_type == VertexType::Internal
                        && (self.graph.edges_directed(*index, Incoming).count() != 1
                            || self.graph.edges_directed(*index, Outgoing).count() != 1)
                })
                .collect();
            let num_bus_nodes = bus_indices.len();

            // The coefficient matrix for the system of equations
            let mut coeffs: Mat<T> = Mat::zeros(num_paths + num_bus_nodes, num_unknowns);
            // The column vector of RHS values for the system of equations
            let mut column: Col<T> = Col::zeros(num_paths + num_bus_nodes);

            // Begin by establishing equations for the voltage drop along every source->sink
            // path. (KVL)
            for (i, path) in paths.iter().enumerate() {
                // These unwraps are safe since we know the paths exist and since source nodes
                // are guaranteed to have a voltage recorded.
                column.write(i, self.graph.node_weight(path.0).unwrap().voltage.unwrap());

                for edge_index in &path.1 {
                    let edge = self.graph.edge_weight(*edge_index).unwrap();
                    let current_index = edge.current_id.unwrap();
                    // Add on this edge's contribution to what's there, it might be in series
                    let current_val = coeffs.read(i, current_index);
                    coeffs.write(i, current_index, edge.admittance.faer_inv() + current_val);
                }
            }

            // For each of the bussing nodes, create an equation from the fact that current
            // cannot pool anywhere, whatever enters a node must also leave. (KCL)
            for (i, node_index) in (num_paths..(num_paths + num_bus_nodes)).zip(bus_indices) {
                // Place a +1 coefficient on each incoming current
                for edge_ref in self.graph.edges_directed(node_index, Incoming) {
                    let current_index = edge_ref.weight().current_id.unwrap();
                    coeffs.write(i, current_index, T::faer_one());
                }
                // Place a -1 coefficient on each outgoing current
                for edge_ref in self.graph.edges_directed(node_index, Outgoing) {
                    let current_index = edge_ref.weight().current_id.unwrap();
                    coeffs.write(i, current_index, T::faer_one().faer_neg());
                }
            }

            // Solve the system of equations
            let solver = Svd::new(coeffs.as_ref());
            solver.solve_in_place(column.as_mut().as_2d_mut());
            column.resize_with(num_unknowns, |_| T::faer_zero());

            for edge_weight in self.graph.edge_weights_mut() {
                let current_id = edge_weight.current_id.unwrap();
                edge_weight.current = Some(column.read(current_id));
            }
        }

        /// Solve for the voltages at every node in the circuit. Uses the currents
        /// assigned to each [`EdgeMetadata`] by `solve_currents`.
        ///
        /// Returns nothing as the voltages are set on each node's [`VertexMetadata`].
        pub fn solve_voltages(&mut self) {
            // First ensure that the currents have in fact been found. If they aren't,
            // it's not possible to find the voltages.
            if self
                .graph
                .edge_weights()
                .map(|weight| weight.current)
                .any(|current| current.is_none())
            {
                self.solve_currents()
            }

            // We ignore the possible error here since for it to occur, there would
            // need to be a cycle in the graph, which would have mucked up finding
            // the currents in the first place, and shouldn't occur anyway.
            let sorted_nodes = algo::toposort(&self.graph, None).unwrap();

            for node_index in sorted_nodes {
                let weight = self.graph.node_weight(node_index).unwrap();

                if weight.vertex_type != VertexType::Internal {
                    continue;
                }

                let prior_index = self
                    .graph
                    .neighbors_directed(node_index, Incoming)
                    .next()
                    .unwrap();
                let prior_voltage = self
                    .graph
                    .node_weight(prior_index)
                    .unwrap()
                    .voltage
                    .unwrap();
                let connection_index = self.graph.find_edge(prior_index, node_index).unwrap();
                let edge_weight = self.graph.edge_weight(connection_index).unwrap();

                let new_voltage = Some(
                    prior_voltage
                        - edge_weight.current.unwrap() * edge_weight.admittance.faer_inv(),
                );

                let weight = self.graph.node_weight_mut(node_index).unwrap();

                weight.voltage = new_voltage;
            }
        }

        /// Find the power consumed by an edge.
        ///
        /// For complex-valued systems, the returned value will be in the form
        /// `P + jQ` where
        /// - `P` is real power,
        /// - `Q` is reactive power, and
        /// - `j` is the imaginary unit.
        ///
        /// Returns [`None`] iff the provided `edge_index` doesn't exist in the graph.
        pub fn power_on_edge(&mut self, edge_index: EdgeIndex) -> Option<T> {
            let raw_edge = self.graph.edge_weight(edge_index);

            if raw_edge.is_none() {
                return None;
            }

            let edge = raw_edge.unwrap();

            if !edge.power.is_none() {
                return edge.power;
            }

            if edge.current.is_none() {
                self.solve_currents();
            }

            let edge = self.graph.edge_weight(edge_index).unwrap();

            let power =
                edge.current.unwrap().faer_mul(edge.current.unwrap()) * edge.admittance.faer_inv();

            let edge = self.graph.edge_weight_mut(edge_index).unwrap();
            edge.power = Some(power);

            Some(power)
        }

        /// Find the power available at any given node.
        ///
        /// For complex-valued systems, the returned value will be in the form
        /// `P + jQ` where
        /// - `P` is real power,
        /// - `Q` is reactive power, and
        /// - `j` is the imaginary unit.
        ///
        /// Returns [`None`] iff the provided `node_index` does not actually exist
        /// in the graph.
        pub fn power_at_node(&mut self, node_index: NodeIndex) -> Option<T> {
            let raw_node = self.graph.node_weight(node_index);

            if raw_node.is_none() {
                return None;
            }

            let node = raw_node.unwrap();

            if !node.power.is_none() {
                return node.power;
            }

            if node.voltage.is_none() {
                self.solve_voltages();
            }

            let node = self.graph.node_weight(node_index).unwrap();

            let power = node.voltage.unwrap()
                * self
                    .graph
                    .edges_directed(node_index, Outgoing)
                    .map(|edge| edge.weight().current.unwrap())
                    .fold(T::faer_zero(), |a, b| a + b);

            let node = self.graph.node_weight_mut(node_index).unwrap();
            node.power = Some(power);

            Some(power)
        }

        /// Eagerly compute the power consumption on every edge, and the power
        /// available at every node.
        ///
        /// The resulting values are stored on the [`VertexMetadata`] and
        /// [`EdgeMetadata`] objects, but can be accessed through the
        /// `current_on_edge()` and `current_at_node()` methods.
        pub fn compute_power(&mut self) {
            for edge_index in self.graph.edge_indices() {
                _ = self.power_on_edge(edge_index);
            }

            for node_index in self.graph.node_indices() {
                _ = self.power_at_node(node_index);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::circuit_graph::*;

    use faer_core::c64;
    use faer_entity::ComplexField;

    #[test]
    fn create_single_vertex() {
        let v: VertexMetadata<f64> = VertexMetadata::new(Some(1f64), 0, VertexType::Source);

        assert!(v.voltage.unwrap() == 1f64);
        assert!(v.tag == 0);
        assert!(v.vertex_type == VertexType::Source);
    }

    #[test]
    fn create_single_edge() {
        let e: EdgeMetadata<f64> = EdgeMetadata::new(0, 1, 0.5);

        assert!(e.tail == 0);
        assert!(e.head == 1);
        assert!(e.admittance == 0.5);
    }

    #[test]
    fn create_graph() {
        let v1 = VertexMetadata::new(Some(1.0), 1, VertexType::Source);
        let v2 = VertexMetadata::new(Some(0.0), 2, VertexType::Sink);

        let e = EdgeMetadata::new(1, 2, 2.0);

        let circuit = Circuit::new(vec![v1, v2], vec![e]);

        assert!(circuit.graph.node_count() == 2);
    }

    /// Set up a simple circuit:
    ///
    /// ```raw
    /// ----------------
    /// |              |
    /// |+            | |
    /// V 3       2ohm|R|
    /// |-            | |
    /// |              |
    /// ----------------
    /// ```
    fn create_simple_circuit() -> Circuit<f64> {
        let source = VertexMetadata::new(Some(3.0), 0, VertexType::Source);
        let sink = VertexMetadata::new(Some(0.0), 1, VertexType::Sink);

        let edge = EdgeMetadata::new(0, 1, 0.5);

        Circuit::new(vec![source, sink], vec![edge])
    }

    /// Test that the circuit was set up properly in that the voltage of the
    /// voltage source was set to the correct value.
    #[test]
    fn check_simple_voltage_source() {
        let circuit = create_simple_circuit();

        let source = circuit
            .graph
            .node_weights()
            .find(|v| v.vertex_type == VertexType::Source)
            .unwrap();
        assert!(source.voltage.unwrap() == 3.0);
    }

    /// Test that the correct number of unknown currents and voltages are being reported.
    #[test]
    fn test_simple_num_unknowns() {
        let mut circuit = create_simple_circuit();

        assert_eq!(circuit.determine_unknown_currents(), 1);
        assert_eq!(circuit.count_unknown_voltages(), 0);
    }

    /// Test that the correct number of paths are found.
    #[test]
    fn test_simple_num_paths() {
        let circuit = create_simple_circuit();

        assert_eq!(circuit.find_paths().len(), 1);
    }

    /// Test that the solved current value through the resistor is correct.
    #[test]
    fn test_simple_solved_current() {
        let mut circuit = create_simple_circuit();

        circuit.solve_currents();
        let solved_currents: Vec<f64> = circuit
            .graph
            .edge_weights()
            .map(|weight| weight.current.unwrap())
            .collect();

        assert_eq!(solved_currents.len(), 1);
        assert!(solved_currents[0] - 1.5 < 1e-10);
    }

    /// Set up a slightly more complex circuit:
    ///
    /// ```raw
    ///     _____
    /// ----__R__---------------
    /// |   1ohm   |           |
    /// |+        | |         | |
    /// V 2   1ohm|R|   0.5ohm|R|
    /// |-        | |         | |
    /// |          |           |
    /// ------------------------
    /// ```
    fn create_complex_circuit() -> Circuit<f64> {
        let source = VertexMetadata::new(Some(2.0), 0, VertexType::Source);

        let v1 = VertexMetadata::new(None, 1, VertexType::Internal);

        let sink = VertexMetadata::new(Some(0.0), 2, VertexType::Sink);

        let e1 = EdgeMetadata::new(0, 1, 1.0);
        let e2 = EdgeMetadata::new(1, 2, 1.0);
        let e3 = EdgeMetadata::new(1, 2, 2.0);

        Circuit::new(vec![source, v1, sink], vec![e1, e2, e3])
    }

    /// Test that the circuit's voltage source was set to the correct value.
    #[test]
    fn check_complex_voltage_source() {
        let circuit = create_complex_circuit();

        let source = circuit
            .graph
            .node_weights()
            .find(|v| v.vertex_type == VertexType::Source)
            .unwrap();

        assert_eq!(source.voltage.unwrap(), 2.0);
    }

    /// Test that the correct number of unknowns are reported.
    #[test]
    fn test_complex_num_unknowns() {
        let mut circuit = create_complex_circuit();

        assert_eq!(circuit.determine_unknown_currents(), 3);
        assert_eq!(circuit.count_unknown_voltages(), 1);
    }

    /// Test that the correct number of paths are found.
    #[test]
    fn test_complex_num_paths() {
        let circuit = create_complex_circuit();

        assert_eq!(circuit.find_paths().len(), 2);
    }

    /// Test that the solved current value through each resistor is correct.
    #[test]
    fn test_complex_solved_currents() {
        let mut circuit = create_complex_circuit();

        circuit.solve_currents();
        let solved_currents: Vec<f64> = circuit
            .graph
            .edge_weights()
            .map(|weight| weight.current.unwrap())
            .collect();

        assert_eq!(solved_currents.len(), 3);
        assert!(solved_currents[0] - 1.5 < 1e-10);
        assert!(solved_currents[1] - 0.5 < 1e-10);
        assert!(solved_currents[2] - 1.0 < 1e-10);
    }

    /// Test that the solved voltage drop across each resistor is correct.
    #[test]
    fn test_complex_solved_voltages() {
        let mut circuit = create_complex_circuit();

        circuit.solve_currents();
        circuit.solve_voltages();

        let node = circuit
            .graph
            .node_weights()
            .find(|node| node.vertex_type == VertexType::Internal)
            .unwrap();

        assert!(node.voltage.unwrap() - 0.5 < 1e-10);
    }

    /// Set up a more complex circuit with series components
    ///
    /// ```raw
    ///     _____       _____
    /// ----__R__-------__R__---
    /// |   1ohm   |    0.5ohm |
    /// |+        | |         | |
    /// V 2   1ohm|R|     2ohm|R|
    /// |-        | |         | |
    /// |          |           |
    /// ------------------------
    /// ```
    fn create_series_circuit() -> Circuit<f64> {
        let source = VertexMetadata::new(Some(2.0), 0, VertexType::Source);

        let v1 = VertexMetadata::new(None, 1, VertexType::Internal);
        let v2 = VertexMetadata::new(None, 2, VertexType::Internal);

        let sink = VertexMetadata::new(Some(0.0), 3, VertexType::Sink);

        let e1 = EdgeMetadata::new(0, 1, 1.0);
        let e2 = EdgeMetadata::new(1, 3, 1.0);
        let e3 = EdgeMetadata::new(1, 2, 2.0);
        let e4 = EdgeMetadata::new(2, 3, 0.5);

        Circuit::new(vec![source, v1, v2, sink], vec![e1, e2, e3, e4])
    }

    /// Test that the correct number of unknowns are being reported.
    #[test]
    fn test_series_num_unknowns() {
        let mut circuit = create_series_circuit();

        assert_eq!(circuit.determine_unknown_currents(), 3);
        assert_eq!(circuit.count_unknown_voltages(), 2);
    }

    /// Test that the correct number of paths are found.
    #[test]
    fn test_series_num_paths() {
        let circuit = create_series_circuit();

        assert_eq!(circuit.find_paths().len(), 2);
    }

    /// Test that the correct current values have been found
    #[test]
    fn test_series_solved_currents() {
        let mut circuit = create_series_circuit();

        circuit.solve_currents();
        let solved_currents: Vec<f64> = circuit
            .graph
            .edge_weights()
            .map(|weight| weight.current.unwrap())
            .collect();

        assert_eq!(solved_currents.len(), 4);

        assert!(solved_currents[0] - 7.0 / 6.0 < 1e-10);
        assert!(solved_currents[1] - 5.0 / 6.0 < 1e-10);
        assert!(solved_currents[2] - 1.0 / 3.0 < 1e-10);
        assert!(solved_currents[3] - 1.0 / 3.0 < 1e-10);
    }

    /// Test that the correct voltage values have been found
    #[test]
    fn test_series_solved_voltages() {
        let mut circuit = create_series_circuit();

        circuit.solve_currents();
        circuit.solve_voltages();

        let voltages: Vec<f64> = circuit
            .graph
            .node_weights()
            .map(|weight| weight.voltage.unwrap())
            .collect();

        assert!(voltages[0] - 2.0 < 1e-10);
        assert!(voltages[1] - 5.0 / 6.0 < 1e-10);
        assert!(voltages[2] - 2.0 / 3.0 < 1e-10);
        assert!(voltages[3] - 0.0 < 1e-10);
    }

    /// Set up a circuit with multiple source and sink nodes. The location of ground
    /// for analysis is marked with G.
    ///
    /// ```raw
    ///        _____
    /// -------__R__--------------------------
    /// |      1ohm        |                 |
    /// |+                | |               | |
    /// |V 3          1ohm|R|           1ohm|R|
    /// |-                | |               | |
    /// |      _____       |                 |
    /// |------__R__-------------------------|
    /// |      1ohm        |                 |
    /// |+                | |               | |
    /// |V 2          1ohm|R|           1ohm|R|
    /// |-                | |               | |
    /// |                  |                 |
    /// --------------------------------------
    ///          |
    ///          G
    /// ```
    fn create_multiple_source_circuit() -> Circuit<f64> {
        // Note that source0 represents the node at the very top left of the circuit,
        // and thus has the combined voltage of both the circuit's voltage sources.
        let source0 = VertexMetadata::new(Some(5.0), 0, VertexType::Source);
        let source1 = VertexMetadata::new(Some(2.0), 1, VertexType::Source);

        let v0 = VertexMetadata::new(None, 2, VertexType::Internal);
        let v1 = VertexMetadata::new(None, 3, VertexType::Internal);

        let sink = VertexMetadata::new(Some(0.0), 4, VertexType::Sink);

        let e0 = EdgeMetadata::new(0, 2, 1.0);
        let e1 = EdgeMetadata::new(2, 3, 1.0);
        let e2 = EdgeMetadata::new(2, 3, 1.0);
        let e3 = EdgeMetadata::new(1, 3, 1.0);
        let e4 = EdgeMetadata::new(3, 4, 1.0);
        let e5 = EdgeMetadata::new(3, 4, 1.0);

        Circuit::new(
            vec![source0, source1, v0, v1, sink],
            vec![e0, e1, e2, e3, e4, e5],
        )
    }

    /// Test that the correct number of paths are found.
    #[test]
    fn test_multiple_num_paths() {
        let circuit = create_multiple_source_circuit();

        assert_eq!(circuit.find_paths().len(), 6);
    }

    /// Test that the correct currents are found.
    #[test]
    fn test_multiple_solved_currents() {
        let mut circuit = create_multiple_source_circuit();

        circuit.solve_currents();
        let solved_currents: Vec<f64> = circuit
            .graph
            .edge_weights()
            .map(|weight| weight.current.unwrap())
            .collect();

        assert_eq!(solved_currents.len(), 6);

        assert!(solved_currents[0] - 26. / 11. < 1e-10);
        assert!(solved_currents[1] - 13. / 11. < 1e-10);
        assert!(solved_currents[2] - 13. / 11. < 1e-10);
        assert!(solved_currents[3] - 6. / 11. < 1e-10);
        assert!(solved_currents[4] - 16. / 11. < 1e-10);
        assert!(solved_currents[5] - 16. / 11. < 1e-10);
    }

    /// Test that the correct voltages are found.
    #[test]
    fn test_multiple_solved_voltages() {
        let mut circuit = create_multiple_source_circuit();

        circuit.solve_currents();
        circuit.solve_voltages();

        let voltages: Vec<f64> = circuit
            .graph
            .node_weights()
            .map(|weight| weight.voltage.unwrap())
            .collect();

        assert!(voltages[0] - 5.0 < 1e-10);
        assert!(voltages[1] - 2.0 < 1e-10);
        assert!(voltages[2] - 29.0 / 11.0 < 1e-10);
        assert!(voltages[3] - 16.0 / 11.0 < 1e-10);
        assert!(voltages[4] - 0.0 < 1e-10);
    }

    /// Set up an AC circuit.
    ///
    /// ```raw
    /// ------|C(------
    /// |   -j2ohm    |
    /// |+            $
    /// |V 4    j1ohm I
    /// |-            $
    /// |    _____    |
    /// -----__R__-----
    ///       2ohm
    /// ```
    fn create_ac_circuit() -> Circuit<c64> {
        let source = VertexMetadata::new(Some(c64::new(4.0, 0.0)), 0, VertexType::Source);
        let v1 = VertexMetadata::new(None, 1, VertexType::Internal);
        let v2 = VertexMetadata::new(None, 2, VertexType::Internal);
        let sink = VertexMetadata::new(Some(c64::faer_zero()), 3, VertexType::Sink);

        let e1 = EdgeMetadata::new(0, 1, c64::new(0.0, 0.5));
        let e2 = EdgeMetadata::new(1, 2, c64::new(0.0, -1.0));
        let e3 = EdgeMetadata::new(2, 3, c64::new(0.5, 0.0));

        Circuit::new(vec![source, v1, v2, sink], vec![e1, e2, e3])
    }

    /// Test that the correct currents are found.
    #[test]
    fn test_ac_currents() {
        let mut circuit = create_ac_circuit();

        circuit.solve_currents();

        let solved_currents: Vec<c64> = circuit
            .graph
            .edge_weights()
            .map(|weight| weight.current.unwrap())
            .collect();

        for current in solved_currents {
            assert!((current - c64::new(8.0 / 5.0, 4.0 / 5.0)).faer_abs() < 1e-10);
        }
    }

    /// Test that the correct voltages are found.
    #[test]
    fn test_ac_voltages() {
        let mut circuit = create_ac_circuit();

        circuit.solve_voltages();

        let solved_voltages: Vec<c64> = circuit
            .graph
            .node_weights()
            .map(|weight| weight.voltage.unwrap())
            .collect();

        assert!((solved_voltages[0] - c64::new(4.0, 0.0)).faer_abs() < 1e-10);
        assert!((solved_voltages[1] - c64::new(12.0 / 5.0, 16.0 / 5.0)).faer_abs() < 1e-10);
        assert!((solved_voltages[2] - c64::new(16.0 / 5.0, 8.0 / 5.0)).faer_abs() < 1e-10);
        assert!((solved_voltages[3] - c64::faer_zero()).faer_abs() < 1e-10);
    }
}

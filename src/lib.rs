pub mod circuit_graph {
    use std::collections::HashMap;

    use faer::prelude::SpSolver;
    use faer::solvers::FullPivLu;
    use faer::Mat;
    use faer_entity::ComplexField;
    use petgraph::prelude::*;

    #[derive(PartialEq, Eq)]
    pub enum VertexType {
        Source,
        Sink,
        Internal,
    }

    pub struct VertexMetadata<T> {
        pub voltage: T,
        pub tag: u32,
        pub vertex_type: VertexType,
    }

    impl<T> VertexMetadata<T> {
        pub fn new(voltage: T, tag: u32, vertex_type: VertexType) -> Self {
            Self {
                voltage,
                tag,
                vertex_type,
            }
        }
    }

    pub struct EdgeMetadata<T> {
        pub tail: u32,
        pub head: u32,
        pub conductance: T,
        current_id: Option<usize>,
    }

    impl<T> EdgeMetadata<T> {
        pub fn new(tail: u32, head: u32, conductance: T) -> Self {
            Self {
                tail,
                head,
                conductance,
                current_id: None,
            }
        }
    }

    /// A struct representing a passive circuit. It relies on a flow graph where
    ///  - vertices are annotated with a [`VertexType`]
    ///  - vertices are annotated with a voltage
    ///  - edges are annotated with a weight equal to the conductance of the
    /// component (reciprocal of resistance).
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
            // resistor whose current we want to find, but that resistors in series (which
            // can be uniquely associated with vertices having deg_in = deg_out = 1)
            // need to be uncounted since any group all share the same current.
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

            println!("{}", next_current_index);

            let out = self.graph.edge_count() - series_nodes.len();

            println!("{}", out);

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

    impl<T: ComplexField + std::ops::Add<Output = T>> Circuit<T> {
        /// Solve current values
        pub fn solve_currents(&mut self) -> Vec<T> {
            let num_unkowns = self.determine_unknown_currents();

            let mut coeffs: Mat<T> = Mat::zeros(num_unkowns, num_unkowns);
            let mut voltages: Mat<T> = Mat::zeros(num_unkowns, 1);

            let paths = self.find_paths();
            let num_paths = paths.len();

            // Begin by establishing equations for the voltage drop along every source->sink
            // path
            for (i, path) in paths.iter().enumerate() {
                voltages.write(i, 0, self.graph.node_weight(path.0).unwrap().voltage);

                for edge_index in &path.1 {
                    let edge = self.graph.edge_weight(*edge_index).unwrap();
                    let current_index = edge.current_id.unwrap();
                    // Add on this edge's contribution to what's there, it might be in series
                    let current_val = coeffs.read(i, current_index);
                    coeffs.write(i, current_index, edge.conductance.faer_inv() + current_val);
                }
            }

            // Now find internal vertices which bus multiple branches
            let indices: Vec<NodeIndex> = self
                .graph
                .node_indices()
                .filter(|index| {
                    self.graph.node_weight(*index).unwrap().vertex_type == VertexType::Internal
                        && (self.graph.edges_directed(*index, Incoming).count() != 1
                            || self.graph.edges_directed(*index, Outgoing).count() != 1)
                })
                .collect();

            // For each of these vertices, create an equation from the fact that current
            // cannot pool anywhere, whatever enters a node must also leave.
            for (i, node_index) in (num_paths..num_unkowns).zip(indices) {
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
            let solver = FullPivLu::new(coeffs.as_ref());
            let result = solver.solve(voltages);

            let mut out = Vec::new();

            for i in 0..num_unkowns {
                out.push(result.read(i, 0));
            }

            out
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::circuit_graph::*;

    #[test]
    fn create_single_vertex() {
        let v: VertexMetadata<f64> = VertexMetadata::new(1f64, 0, VertexType::Source);

        assert!(v.voltage == 1f64);
        assert!(v.tag == 0);
        assert!(v.vertex_type == VertexType::Source);
    }

    #[test]
    fn create_single_edge() {
        let e: EdgeMetadata<f64> = EdgeMetadata::new(0, 1, 0.5);

        assert!(e.tail == 0);
        assert!(e.head == 1);
        assert!(e.conductance == 0.5);
    }

    #[test]
    fn create_graph() {
        let v1 = VertexMetadata::new(1.0, 1, VertexType::Source);
        let v2 = VertexMetadata::new(0.0, 2, VertexType::Sink);

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
    /// V             |R|
    /// |-            | |
    /// |              |
    /// ----------------
    /// ```
    fn create_simple_circuit() -> Circuit<f64> {
        let source = VertexMetadata::new(3.0, 0, VertexType::Source);
        let sink = VertexMetadata::new(0.0, 1, VertexType::Sink);

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
        assert!(source.voltage == 3.0);
    }

    /// Test that the correct number of unknown currents and voltages are being reported.
    #[test]
    fn check_simple_num_unknowns() {
        let mut circuit = create_simple_circuit();

        assert_eq!(circuit.determine_unknown_currents(), 1);
        assert_eq!(circuit.count_unknown_voltages(), 0);
    }

    /// Test that the correct number of paths are found.
    #[test]
    fn check_simple_num_paths() {
        let circuit = create_simple_circuit();

        assert_eq!(circuit.find_paths().len(), 1);
    }

    /// Test that the solved current value through the resistor is correct.
    #[test]
    fn test_simple_solved_current() {
        let mut circuit = create_simple_circuit();

        let solved_currents = circuit.solve_currents();

        assert_eq!(solved_currents.len(), 1);
        assert!(solved_currents[0] - 1.5 < 1e-10);
    }

    /// Test that the solved voltage drop across the resistor is correct.
    #[test]
    fn test_simple_solved_voltage() {
        todo!()
    }

    /// Set up a slightly more complex circuit:
    ///
    /// ```raw
    ///     __ __
    /// ----__R__---------------
    /// |   1ohm   |           |
    /// |+        | |         | |
    /// V 2   1ohm|R|   0.5ohm|R|
    /// |-        | |         | |
    /// |          |           |
    /// ------------------------
    /// ```
    fn create_complex_circuit() -> Circuit<f64> {
        let source = VertexMetadata::new(2.0, 0, VertexType::Source);

        let v1 = VertexMetadata::new(-1.0, 1, VertexType::Internal);

        let sink = VertexMetadata::new(0.0, 2, VertexType::Sink);

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

        assert_eq!(source.voltage, 2.0);
    }

    /// Test that the correct number of unknowns are reported.
    #[test]
    fn check_complex_num_unknowns() {
        let mut circuit = create_complex_circuit();

        assert_eq!(circuit.determine_unknown_currents(), 3);
        assert_eq!(circuit.count_unknown_voltages(), 1);
    }

    /// Test that the correct number of paths are found.
    #[test]
    fn check_complex_num_paths() {
        let circuit = create_complex_circuit();

        assert_eq!(circuit.find_paths().len(), 2);
    }

    /// Test that the solved current value through each resistor is correct.
    #[test]
    fn test_complex_solved_currents() {
        let mut circuit = create_complex_circuit();

        let solved_currents = circuit.solve_currents();

        assert_eq!(solved_currents.len(), 3);
        assert!(solved_currents[0] - 1.5 < 1e-10);
        assert!(solved_currents[1] - 0.5 < 1e-10);
        assert!(solved_currents[2] - 1.0 < 1e-10);
    }

    /// Test that the solved voltage drop across each resistor is correct.
    #[test]
    fn test_complex_solved_voltages() {
        todo!()
    }

    /// Set up a more complex circuit with series components
    ///
    /// ```raw
    ///     __ __       __ __
    /// ----__R__-------__R__---
    /// |   1ohm   |    0.5ohm |
    /// |+        | |         | |
    /// V     1ohm|R|     2ohm|R|
    /// |-        | |         | |
    /// |          |           |
    /// ------------------------
    /// ```
    fn create_series_circuit() -> Circuit<f64> {
        let source = VertexMetadata::new(2.0, 0, VertexType::Source);

        let v1 = VertexMetadata::new(-1.0, 1, VertexType::Internal);
        let v2 = VertexMetadata::new(-1.0, 2, VertexType::Internal);

        let sink = VertexMetadata::new(0.0, 3, VertexType::Sink);

        let e1 = EdgeMetadata::new(0, 1, 1.0);
        let e2 = EdgeMetadata::new(1, 3, 1.0);
        let e3 = EdgeMetadata::new(1, 2, 2.0);
        let e4 = EdgeMetadata::new(2, 3, 0.5);

        Circuit::new(vec![source, v1, v2, sink], vec![e1, e2, e3, e4])
    }

    /// Test that the correct number of unknowns are being reported.
    #[test]
    fn check_series_num_unknowns() {
        let mut circuit = create_series_circuit();

        assert_eq!(circuit.determine_unknown_currents(), 3);
        assert_eq!(circuit.count_unknown_voltages(), 2);
    }

    /// Test that the correct number of paths are found.
    #[test]
    fn check_series_num_paths() {
        let circuit = create_series_circuit();

        assert_eq!(circuit.find_paths().len(), 2);
    }

    /// Test that the correct current values have been found
    #[test]
    fn check_series_solved_currents() {
        let mut circuit = create_series_circuit();

        let solved_currents = circuit.solve_currents();

        assert_eq!(solved_currents.len(), 3);

        assert!(solved_currents[0] - 1.0 / 3.0 < 1e-10);
        assert!(solved_currents[1] - 7.0 / 6.0 < 1e-10);
        assert!(solved_currents[2] - 5.0 / 6.0 < 1e-10);
    }

    /// Set up a circuit with multiple source and sink nodes.
    ///
    /// ```raw
    ///
    /// -------__R__--------------__R__-------
    /// |                  |                 |
    /// |                 | |                |
    /// |+                |R|               | |
    /// |V                | |               |R|
    /// |-                 |                | |
    /// |     ----__R__---------__R__---     |
    /// |     |            |           |     |
    /// |     |+           |           |     |
    /// |     |V           |          | |    |
    /// |     |-           |          |R|    |
    /// |     |            |          | |    |
    /// |     ----__R__-----           |     |
    /// |                              |     |
    /// --------------------------------------
    /// ```
    fn create_multiple_source_circuit() -> Circuit<f64> {
        let source1 = VertexMetadata::new(2.0, 0, VertexType::Source);
        let source2 = VertexMetadata::new(1.5, 1, VertexType::Source);

        let v1 = VertexMetadata::new(-1.0, 2, VertexType::Internal);
        let v2 = VertexMetadata::new(-1.0, 3, VertexType::Internal);
        let v3 = VertexMetadata::new(-1.0, 4, VertexType::Internal);
        let v4 = VertexMetadata::new(-1.0, 5, VertexType::Internal);

        let sink1 = VertexMetadata::new(0.0, 6, VertexType::Sink);
        let sink2 = VertexMetadata::new(0.0, 7, VertexType::Sink);

        let e1 = EdgeMetadata::new(0, 2, 1.0);
        let e2 = EdgeMetadata::new(2, 3, 1.0);
        let e3 = EdgeMetadata::new(3, 6, 1.0);
        let e4 = EdgeMetadata::new(2, 4, 1.0);
        let e5 = EdgeMetadata::new(4, 5, 1.0);
        let e6 = EdgeMetadata::new(5, 6, 1.0);
        let e7 = EdgeMetadata::new(4, 7, 1.0);
        let e8 = EdgeMetadata::new(1, 4, 1.0);

        Circuit::new(
            vec![source1, source2, v1, v2, v3, v4, sink1, sink2],
            vec![e1, e2, e3, e4, e5, e6, e7, e8],
        )
    }

    /// Test that the correct number of paths are found.
    #[test]
    fn check_multiple_num_paths() {
        let circuit = create_multiple_source_circuit();

        assert_eq!(circuit.find_paths().len(), 5);
    }
}

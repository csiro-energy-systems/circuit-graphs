pub mod circuit_graph {
    use std::collections::{HashMap, HashSet};

    use faer::solvers::{FullPivLu, SolverCore};
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
    }

    impl<T> EdgeMetadata<T> {
        pub fn new(tail: u32, head: u32, conductance: T) -> Self {
            Self {
                tail,
                head,
                conductance,
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
                let tag = vertex.tag.clone();
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

        /// Count the number of unknown currents that need to be found when solving.
        pub fn count_unknown_currents(&self) -> usize {
            // This method uses the fact that every edge of the graph corresponds to a
            // resistor whose current we want to find, but that resistors in series (which
            // can be uniquely associated with vertices having deg_in = deg_out = 1)
            // need to be uncounted since any group all share the same current.
            self.graph.edge_count()
                - self
                    .graph
                    .node_indices()
                    .filter(|index| {
                        self.graph.edges_directed(*index, Incoming).count() == 1
                            && self.graph.edges_directed(*index, Outgoing).count() == 1
                    })
                    .count()
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

        /// Finds, as sequences of edges, all the paths between a given source and sink node.
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

    impl<T: ComplexField> Circuit<T> {
        /// Solve current values
        pub fn solve_currents(&self) {
            let num_unkowns = self.count_unknown_currents();

            let mut coeffs: Mat<T> = Mat::zeros(num_unkowns, num_unkowns);
            let mut voltages: Mat<T> = Mat::zeros(num_unkowns, 1);

            let paths = self.find_paths();
            let num_paths = paths.len();

            // Begin by establishing equations for the current drop along every source->sink
            // path
            for (i, path) in paths.iter().enumerate() {
                voltages.write(i, 0, self.graph.node_weight(path.0).unwrap().voltage);

                for edge_index in &path.1 {
                    let edge = self.graph.edge_weight(*edge_index).unwrap();
                    coeffs.write(i, edge_index.index(), edge.conductance.faer_inv());
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
                    coeffs.write(i, edge_ref.id().index(), T::faer_one());
                }
                // Place a -1 coefficient on each outgoing current
                for edge_ref in self.graph.edges_directed(node_index, Outgoing) {
                    coeffs.write(i, edge_ref.id().index(), T::faer_one().faer_neg());
                }
            }

            println!("coeffs: \n {:#?} \n voltages: \n {:#?}", coeffs, voltages);

            let solver = FullPivLu::new(coeffs.as_ref());
            let inverse = solver.inverse();

            let result = inverse * voltages;

            println!("result: {:#?}", result);
        }
    }

    /// A struct to hold the functionality for determining sets of edges in series
    /// which will share current. This way we can figure out a mapping from edge
    /// index to current value.
    struct SeriesSets {
        pub sets: Vec<HashSet<EdgeIndex>>,
    }

    impl SeriesSets {
        pub fn new<T>(circuit: &Circuit<T>) -> Self {
            let mut new = Self { sets: Vec::new() };

            // Essentially don't bother with the hard bit if we don't need to
            if circuit.count_unknown_currents() == circuit.graph.edge_count() {
                for edge_index in circuit.graph.edge_indices() {
                    new.insert(edge_index);
                }

                return new;
            }

            let indices = circuit.graph.node_indices().filter(|index| {
                circuit.graph.edges_directed(*index, Incoming).count() == 1
                    && circuit.graph.edges_directed(*index, Outgoing).count() == 1
            });

            for node_index in indices {
                // We know because of that filter above that there will be exactly one
                // result in the iterator, so the unwrap is safe.
                let incoming_edge = circuit
                    .graph
                    .edges_directed(node_index, Incoming)
                    .next()
                    .unwrap()
                    .id();
                let outgoing_edge = circuit
                    .graph
                    .edges_directed(node_index, Outgoing)
                    .next()
                    .unwrap()
                    .id();
                new.insert_series_pair(incoming_edge, outgoing_edge);
            }

            for edge_index in circuit.graph.edge_indices() {
                new.insert(edge_index);
            }

            new
        }

        fn insert(&mut self, edge_index: EdgeIndex) {
            // If the edge is already somewhere in here, do nothing.
            if self.sets.iter().any(|set| set.contains(&edge_index)) {
                return;
            }

            // Otherwise, insert in a new HashSet.
            let mut new_set = HashSet::new();
            new_set.insert(edge_index);
            self.sets.push(new_set);
        }

        fn insert_series_pair(&mut self, edge: EdgeIndex, other_edge: EdgeIndex) {
            // .find() is okay here since there will only ever be exactly one HashSet that
            // contains any given edge.
            if let Some(set) = self
                .sets
                .iter_mut()
                .find(|set| set.contains(&edge) || set.contains(&other_edge))
            {
                set.insert(edge);
                set.insert(other_edge);
            } else {
                let mut new_set = HashSet::new();
                new_set.insert(edge);
                new_set.insert(other_edge);
                self.sets.push(new_set);
            }
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
    fn create_simple_circuit(source_voltage: f64, resistance: f64) -> Circuit<f64> {
        let source = VertexMetadata::new(source_voltage, 0, VertexType::Source);
        let sink = VertexMetadata::new(0.0, 1, VertexType::Sink);

        let edge = EdgeMetadata::new(0, 1, 1.0 / resistance);

        Circuit::new(vec![source, sink], vec![edge])
    }

    /// Test that the circuit was set up properly in that the voltage of the
    /// voltage source was set to the correct value.
    #[test]
    fn check_simple_voltage_source() {
        let circuit = create_simple_circuit(5.0, 4.0);

        let source = circuit
            .graph
            .node_weights()
            .find(|v| v.vertex_type == VertexType::Source)
            .unwrap();
        assert!(source.voltage == 5.0);
    }

    /// Test that the correct number of unknown currents and voltages are being reported.
    #[test]
    fn check_simple_num_unknowns() {
        let circuit = create_simple_circuit(2.0, 8.0);

        assert_eq!(circuit.count_unknown_currents(), 1);
        assert_eq!(circuit.count_unknown_voltages(), 0);
    }

    /// Test that the correct number of paths are found.
    #[test]
    fn check_simple_num_paths() {
        let circuit = create_simple_circuit(3.0, 2.0);

        assert_eq!(circuit.find_paths().len(), 1);
    }

    /// Test that the resistance of the resistor was set to the correct value.
    #[test]
    fn check_simple_resistance() {
        todo!()
    }

    /// Test that the solved current value through the resistor is correct.
    #[test]
    fn test_simple_solved_current() {
        todo!()
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
    /// |          |           |
    /// |+        | |         | |
    /// V         |R|         |R|
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
        let circuit = create_complex_circuit();

        assert_eq!(circuit.count_unknown_currents(), 3);
        assert_eq!(circuit.count_unknown_voltages(), 1);
    }

    /// Test that the correct number of paths are found.
    #[test]
    fn check_complex_num_paths() {
        let circuit = create_complex_circuit();

        assert_eq!(circuit.find_paths().len(), 2);
    }

    /// Test that the resistances of the circuit's resistors were set correctly.
    #[test]
    fn check_complex_resistances() {
        todo!()
    }

    /// Test that the solved current value through each resistor is correct.
    #[test]
    fn test_complex_solved_currents() {
        todo!()
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
    /// |          |           |
    /// |+        | |         | |
    /// V         |R|         |R|
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
        let circuit = create_series_circuit();

        assert_eq!(circuit.count_unknown_currents(), 3);
        assert_eq!(circuit.count_unknown_voltages(), 2);
    }

    /// Test that the correct number of paths are found.
    #[test]
    fn check_series_num_paths() {
        let circuit = create_series_circuit();

        assert_eq!(circuit.find_paths().len(), 2);
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

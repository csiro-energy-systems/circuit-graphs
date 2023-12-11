mod circuit_graph {
    use std::collections::HashMap;

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
            self.graph.edge_count()
        }

        /// Count the number of unknown voltages that need to be found when solving.
        pub fn count_unknown_voltages(&self) -> usize {
            self
                .graph
                .node_weights()
                .filter(|v| v.vertex_type == VertexType::Internal)
                .count()
        }
    }
}

#[cfg(test)]
mod test {
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

    /// Setup a simple circuit:
    ///
    /// ```raw
    /// ----------------
    /// |              |
    /// |             | |
    /// V             |R|
    /// |             | |
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

    /// Setup a slightly more complex circuit:
    ///
    /// ```raw
    ///     __ __
    /// ----__R__---------------
    /// |          |           |
    /// |         | |         | |
    /// V         |R|         |R|
    /// |         | |         | |
    /// |          |           |
    /// ------------------------
    /// ```
    fn create_complex_circuit() {
        todo!()
    }

    /// Test that the circuit's voltage source was set to the correct value.
    #[test]
    fn check_complex_voltage_source() {
        todo!()
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
}

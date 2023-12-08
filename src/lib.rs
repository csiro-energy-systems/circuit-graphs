// /// A struct representing a passive circuit. It relies on a flow graph where
// ///  - source vertices are annotated with a voltage
// ///  - edges are annotated with a weight equal to the conductance of the
// /// component (reciprocal of resistance).
// struct Circuit {
//     // Graph doesn't exist, it's just a placeholder for now
//     graph: Graph,
// }

// impl Circuit {
//     /// Create a new circuit
//     fn new() {}

//     /// Solve for the current through and voltage across each of the circuit's
//     /// resistors.
//     fn solve() {}

//     /// Get the source nodes from the graph
//     fn get_sources() {}
// }

mod circuit_graph {
    use std::collections::HashMap;

    use petgraph::{prelude::*, data::Build, graph::Node};

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
            Self { tail, head, conductance }
        }
    }

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
                    edge
                );
            }


            Circuit { graph }
        }
    }
}

#[cfg(test)]
mod test {
    use petgraph::visit::IntoNodeIdentifiers;

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
    /// ----------------
    /// |              |
    /// |             | |
    /// V ^           |R|
    /// |             | |
    /// |              |
    /// ----------------
    fn create_simple_circuit() {
        todo!()
    }

    /// Test that the circuit was set up properly in that the voltage of the
    /// voltage source was set to the correct value.
    #[test]
    fn check_simple_voltage_source() {
        todo!()
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
    ///     __ __
    /// ----__R__---------------
    /// |          |           |
    /// |         | |         | |
    /// V         |R|         |R|
    /// |         | |         | |
    /// |          |           |
    /// ------------------------
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
    fn complex_solved_voltages() {
        todo!()
    }
}

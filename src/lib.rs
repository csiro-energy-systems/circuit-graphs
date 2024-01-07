pub mod circuit_graph {
    use std::collections::HashMap;

    use faer::prelude::SpSolver;
    use faer::solvers::Svd;
    use faer::{Col, Mat};
    use faer_entity::ComplexField;
    use petgraph::algo;
    use petgraph::prelude::*;

    #[derive(PartialEq, Eq)]
    pub enum VertexTag {
        Source {
            tag: u32,
        },
        Sink {
            tag: u32,
        },
        Internal {
            tag: u32,
        },
        Transformer {
            primary_tag: u32,
            primary_coils: u32,
            secondary_tag: u32,
            secondary_coils: u32,
        },
        TransformerPrimary {
            tag: u32,
            secondary_tag: u32,
            num_coils: u32,
        },
        TransformerSecondary {
            tag: u32,
            primary_tag: u32,
            num_coils: u32,
        },
    }

    impl VertexTag {
        /// Returns the tag that this [`VertexTag`] holds, no matter the variant.
        ///
        /// # Panics
        /// This will panic if the [`VertexTag`] is a [`Transformer`][transformer]
        /// variant.
        ///
        /// [transformer]: `VertexTag::Transformer`
        pub fn get_tag(&self) -> &u32 {
            match self {
                Self::Transformer { .. } => {
                    panic!("Attempted to get tag from Transformer node.");
                }
                Self::Source { tag }
                | Self::Sink { tag }
                | Self::Internal { tag }
                | Self::TransformerPrimary {
                    tag,
                    secondary_tag: _,
                    num_coils: _,
                }
                | Self::TransformerSecondary {
                    tag,
                    primary_tag: _,
                    num_coils: _,
                } => tag,
            }
        }

        /// Returns whether this [`VertexTag`] is of the [`Internal`][internal]
        /// variant.
        ///
        /// [internal]: `VertexTag::Internal`
        pub fn is_internal(&self) -> bool {
            match self {
                Self::Internal { .. } => true,
                _ => false,
            }
        }

        /// Returns whether this [`VertexTag`] is of the [`Source`][source]
        /// variant.
        ///
        /// [source]: `VertexTag::Source`
        pub fn is_source(&self) -> bool {
            match self {
                Self::Source { .. } => true,
                _ => false,
            }
        }
    }

    pub struct VertexMetadata<T> {
        pub voltage: Option<T>,
        pub vertex_tag: VertexTag,
        pub power: Option<T>,
    }

    impl<T> VertexMetadata<T> {
        /// Construct a new [`VertexMetadata`]. This constructor guarantees that a
        /// source or sink node will have a voltage associated, and that internal nodes
        /// will have [`None`] instead.
        ///
        /// # Panics
        /// This constructor panics if the given `vertex_tag` is either
        /// `TransformerPrimary` or `TransformerSecondary`; these types are for
        /// internal use and should not be constructed directly.
        pub fn new(voltage: Option<T>, vertex_tag: VertexTag) -> Self {
            match vertex_tag {
                VertexTag::TransformerPrimary { .. } | VertexTag::TransformerSecondary { .. } => {
                    panic!("Attempted to create an internal transformer node via new()");
                }
                VertexTag::Internal { .. } | VertexTag::Transformer { .. } => Self {
                    voltage: None,
                    vertex_tag,
                    power: None,
                },
                _ => {
                    assert!(!voltage.is_none());
                    Self {
                        voltage,
                        vertex_tag,
                        power: None,
                    }
                }
            }
        }
    }

    enum EdgeType<T> {
        Component { admittance: T },
        Transformer,
    }

    /// A struct to hold the information for an edge in the graph.
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
        current_id: Option<usize>,
        pub current: Option<T>,
        pub power: Option<T>,
        edge_type: EdgeType<T>,
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
                current_id: None,
                current: None,
                power: None,
                edge_type: EdgeType::Component { admittance },
            }
        }

        /// Finds the admittance of an edge.
        ///
        /// # Panics
        /// Will panic if called on a [`Transformer`][transformer] edge.
        ///
        /// [transformer]: `EdgeType::Transformer`
        pub fn get_admittance(&self) -> &T {
            match &self.edge_type {
                EdgeType::Component { admittance } => admittance,
                EdgeType::Transformer => {
                    panic!("Attempted to get admittance of a transformer node");
                }
            }
        }

        /// Returns whether this edge is a [`Transformer`][transformer] edge.
        ///
        /// [transformer]: `EdgeType::Transformer`
        pub fn is_transformer(&self) -> bool {
            match self.edge_type {
                EdgeType::Transformer => true,
                EdgeType::Component { admittance: _ } => false,
            }
        }
    }

    /// A struct representing a passive circuit. It relies on a flow graph where
    ///  - vertices are annotated with a [`VertexTag`],
    ///  - vertices are annotated with a voltage (which may be [`None`]), and
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
        /// Create a new circuit from a collection of vertices and edges.
        ///
        /// # Panics
        /// This constructor will panic if an edge's `tail` or `head` refers to a
        /// non-existant node.
        pub fn new(vertices: Vec<VertexMetadata<T>>, mut edges: Vec<EdgeMetadata<T>>) -> Self {
            let mut graph: DiGraph<VertexMetadata<T>, EdgeMetadata<T>> = DiGraph::new();

            let mut vertex_indices: HashMap<u32, NodeIndex> = HashMap::new();

            for vertex in vertices {
                match vertex.vertex_tag {
                    VertexTag::Transformer {
                        primary_tag,
                        primary_coils,
                        secondary_tag,
                        secondary_coils,
                    } => {
                        // For a transformer, we construct the inner nodes and add them.
                        let primary_vertex = VertexMetadata {
                            voltage: None,
                            vertex_tag: VertexTag::TransformerPrimary {
                                tag: primary_tag,
                                secondary_tag,
                                num_coils: primary_coils,
                            },
                            power: None,
                        };
                        let secondary_vertex = VertexMetadata {
                            voltage: None,
                            vertex_tag: VertexTag::TransformerSecondary {
                                tag: secondary_tag,
                                primary_tag,
                                num_coils: secondary_coils,
                            },
                            power: None,
                        };
                        let primary_index = graph.add_node(primary_vertex);
                        let secondary_index = graph.add_node(secondary_vertex);

                        vertex_indices.insert(primary_tag, primary_index);
                        vertex_indices.insert(secondary_tag, secondary_index);

                        // Finally we must also add on the edge we need between the nodes.
                        edges.push(EdgeMetadata {
                            tail: primary_tag,
                            head: secondary_tag,
                            current_id: None,
                            current: None,
                            power: None,
                            edge_type: EdgeType::Transformer,
                        });
                    }
                    _ => {
                        // Any other type of vertex can simply be added into the graph.
                        let tag = *vertex.vertex_tag.get_tag();
                        let node_index = graph.add_node(vertex);
                        vertex_indices.insert(tag, node_index);
                    }
                }
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
                .filter(|v| v.vertex_tag.is_internal())
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
                    self.graph
                        .node_weight(*index)
                        .unwrap()
                        .vertex_tag
                        .is_internal()
                        && (self.graph.edges_directed(*index, Incoming).count() != 1
                            || self.graph.edges_directed(*index, Outgoing).count() != 1)
                })
                .collect();
            let num_bus_nodes = bus_indices.len();

            // Find all the transformer edges
            let transformer_edges: Vec<EdgeIndex> = self
                .graph
                .edge_indices()
                .filter(|index| self.graph.edge_weight(*index).unwrap().is_transformer())
                .collect();
            let num_transformer_edges = transformer_edges.len();

            // The coefficient matrix for the system of equations
            let mut coeffs: Mat<T> = Mat::zeros(
                num_paths + num_bus_nodes + num_transformer_edges,
                num_unknowns,
            );
            // The column vector of RHS values for the system of equations
            let mut column: Col<T> = Col::zeros(num_paths + num_bus_nodes + num_transformer_edges);

            // Begin by establishing equations for the voltage drop along every source->sink
            // path. (KVL)
            for (i, path) in paths.iter().enumerate() {
                // These unwraps are safe since we know the paths exist and since source nodes
                // are guaranteed to have a voltage recorded.
                column.write(i, self.graph.node_weight(path.0).unwrap().voltage.unwrap());

                for edge_index in &path.1 {
                    let edge = self.graph.edge_weight(*edge_index).unwrap();
                    match edge.edge_type {
                        EdgeType::Component { admittance } => {
                            let current_index = edge.current_id.unwrap();
                            // Add on this edge's contribution to what's there, it might be in series
                            let present_value = coeffs.read(i, current_index);
                            coeffs.write(i, current_index, admittance.faer_inv() + present_value);
                        }
                        EdgeType::Transformer => {}
                    }
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

            // At every transformer, we find the current through it by adding up the current
            // coming in and multiplying by the ratio of the coil windings.
            for (i, edge_index) in ((num_paths + num_bus_nodes)
                ..(num_paths + num_bus_nodes + num_transformer_edges))
                .zip(transformer_edges)
            {
                for edge in self
                    .graph
                    .edges_directed(self.graph.edge_endpoints(edge_index).unwrap().0, Incoming)
                {
                    coeffs.write(i, edge.weight().current_id.unwrap(), T::faer_one());
                }
                let (primary_index, secondary_index) =
                    self.graph.edge_endpoints(edge_index).unwrap();
                let primary_coils = match self.graph.node_weight(primary_index).unwrap().vertex_tag
                {
                    VertexTag::TransformerPrimary {
                        tag: _,
                        secondary_tag: _,
                        num_coils,
                    } => num_coils,
                    _ => panic!("Attempted to get primary coils from wrong vertex"),
                };
                let secondary_coils =
                    match self.graph.node_weight(secondary_index).unwrap().vertex_tag {
                        VertexTag::TransformerSecondary {
                            tag: _,
                            primary_tag: _,
                            num_coils,
                        } => num_coils,
                        _ => panic!("Attempted to get secondary coils from wrong vertex"),
                    };
                coeffs.write(
                    i,
                    self.graph
                        .edge_weight(edge_index)
                        .unwrap()
                        .current_id
                        .unwrap(),
                    T::faer_from_f64(f64::from(primary_coils) / f64::from(secondary_coils))
                        .faer_neg(),
                )
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

                // We only need to calculate voltages for internal and transformer nodes.
                match weight.vertex_tag {
                    VertexTag::Internal { tag: _ }
                    | VertexTag::TransformerPrimary {
                        tag: _,
                        secondary_tag: _,
                        num_coils: _,
                    } => {
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
                        let connection_index =
                            self.graph.find_edge(prior_index, node_index).unwrap();
                        let edge_weight = self.graph.edge_weight(connection_index).unwrap();

                        let new_voltage = Some(
                            prior_voltage
                                - edge_weight.current.unwrap()
                                    * edge_weight.get_admittance().faer_inv(),
                        );

                        let weight = self.graph.node_weight_mut(node_index).unwrap();

                        weight.voltage = new_voltage;
                    }
                    // A TransformerSecondary vertex only depends on the primary vertex's voltage
                    // and the coil ratio.
                    VertexTag::TransformerSecondary {
                        tag: _,
                        primary_tag,
                        num_coils: secondary_coils,
                    } => {
                        let primary_weight = self
                            .graph
                            .node_weights()
                            .find(|weight| *weight.vertex_tag.get_tag() == primary_tag)
                            .unwrap();

                        // Ensure we got a TransformerPrimary, if not we can't find the number of
                        // coils.
                        match primary_weight.vertex_tag {
                            VertexTag::TransformerPrimary {
                                tag: _,
                                secondary_tag: _,
                                num_coils: primary_coils,
                            } => {
                                let voltage =
                                    primary_weight.voltage.unwrap().faer_mul(T::faer_from_f64(
                                        f64::from(secondary_coils) / f64::from(primary_coils),
                                    ));

                                let weight = self.graph.node_weight_mut(node_index).unwrap();
                                weight.voltage = Some(voltage);
                            }
                            _ => {
                                panic!("Attempted to find primary transformer vertex, but the VertexType was wrong. Tag: {:?}", primary_tag)
                            }
                        }
                    }
                    _ => {}
                }
            }
        }

        /// Find the power consumed by an edge.
        ///
        /// For complex-valued systems, the returned value will be in the form
        /// `P + jQ` where
        /// - `P` is real power,
        /// - `Q` is reactive power, and
        /// - `j` is the imaginary unit.
        fn find_power_on_edge(&mut self, edge_index: EdgeIndex) {
            let Some(edge) = self.graph.edge_weight(edge_index) else {
                panic!("Tried to access non-existant edge: {:?}", edge_index);
            };

            if edge.current.is_none() {
                self.solve_currents();
            }

            // self was needed as a mutable borrow so we reborrow here
            let edge = self.graph.edge_weight(edge_index).unwrap();

            let power = match edge.edge_type {
                EdgeType::Component { admittance } => {
                    edge.current.unwrap().faer_mul(edge.current.unwrap()) * admittance.faer_inv()
                }
                EdgeType::Transformer => T::faer_zero(),
            };

            let edge = self.graph.edge_weight_mut(edge_index).unwrap();
            edge.power = Some(power);
        }

        /// Find the power available at any given node and store it on the
        /// [`VertexMetadata`]
        ///
        /// For complex-valued systems, the power value will be in the form
        /// `P + jQ` where
        /// - `P` is real power,
        /// - `Q` is reactive power, and
        /// - `j` is the imaginary unit.
        fn find_power_at_node(&mut self, node_index: NodeIndex) {
            let Some(node) = self.graph.node_weight(node_index) else {
                panic!("Tried to access non-existant node: {:?}", node_index);
            };

            if node.voltage.is_none() {
                self.solve_voltages();
            }

            // self was needed as a mutable borrow so we reborrow here
            let node = self.graph.node_weight(node_index).unwrap();

            let power = node.voltage.unwrap()
                * self
                    .graph
                    .edges_directed(node_index, Outgoing)
                    .map(|edge| edge.weight().current.unwrap())
                    .fold(T::faer_zero(), |a, b| a + b);

            let node = self.graph.node_weight_mut(node_index).unwrap();
            node.power = Some(power);
        }

        /// Eagerly compute the power consumption on every edge, and the power
        /// available at every node.
        ///
        /// The resulting values are stored on the [`VertexMetadata`] and
        /// [`EdgeMetadata`] objects.
        pub fn compute_power(&mut self) {
            for edge_index in self.graph.edge_indices() {
                self.find_power_on_edge(edge_index);
            }

            for node_index in self.graph.node_indices() {
                self.find_power_at_node(node_index);
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
        let v: VertexMetadata<f64> = VertexMetadata::new(Some(1f64), VertexTag::Source { tag: 0 });

        assert!(v.voltage.unwrap() == 1f64);
        assert!(v.vertex_tag == VertexTag::Source { tag: 0 });
    }

    #[test]
    fn create_single_edge() {
        let e: EdgeMetadata<f64> = EdgeMetadata::new(0, 1, 0.5);

        assert!(e.tail == 0);
        assert!(e.head == 1);
        assert!(*e.get_admittance() == 0.5);
    }

    #[test]
    fn create_graph() {
        let v1 = VertexMetadata::new(Some(1.0), VertexTag::Source { tag: 1 });
        let v2 = VertexMetadata::new(Some(0.0), VertexTag::Sink { tag: 2 });

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
        let source = VertexMetadata::new(Some(3.0), VertexTag::Source { tag: 0 });
        let sink = VertexMetadata::new(Some(0.0), VertexTag::Sink { tag: 1 });

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
            .find(|v| v.vertex_tag.is_source())
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

    /// Test that the correct power value is found.
    #[test]
    fn test_simple_solved_power() {
        let mut circuit = create_simple_circuit();

        circuit.compute_power();
        let solved_power: Vec<f64> = circuit
            .graph
            .edge_weights()
            .map(|weight| weight.power.unwrap())
            .collect();

        assert!(solved_power[0] - 4.5 < 1e-10);
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
        let source = VertexMetadata::new(Some(2.0), VertexTag::Source { tag: 0 });

        let v1 = VertexMetadata::new(None, VertexTag::Internal { tag: 1 });

        let sink = VertexMetadata::new(Some(0.0), VertexTag::Sink { tag: 2 });

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
            .find(|v| v.vertex_tag.is_source())
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
            .find(|node| node.vertex_tag.is_internal())
            .unwrap();

        assert!(node.voltage.unwrap() - 0.5 < 1e-10);
    }

    /// Test that the correct power drop along each edge and supply at each node
    /// are found.
    #[test]
    fn test_complex_solved_power() {
        let mut circuit = create_complex_circuit();

        circuit.compute_power();

        let edge_powers: Vec<f64> = circuit
            .graph
            .edge_weights()
            .map(|weight| weight.power.unwrap())
            .collect();
        let node_power = circuit
            .graph
            .node_weights()
            .find(|node| node.vertex_tag.is_internal())
            .unwrap()
            .power
            .unwrap();

        assert!(edge_powers[0] - 2.25 < 1e-10);
        assert!(edge_powers[1] - 0.25 < 1e-10);
        assert!(edge_powers[2] - 0.5 < 1e-10);

        assert!(node_power - 0.75 < 1e-10);
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
        let source = VertexMetadata::new(Some(2.0), VertexTag::Source { tag: 0 });

        let v1 = VertexMetadata::new(None, VertexTag::Internal { tag: 1 });
        let v2 = VertexMetadata::new(None, VertexTag::Internal { tag: 2 });

        let sink = VertexMetadata::new(Some(0.0), VertexTag::Sink { tag: 3 });

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

    /// Test that the correct power values have been found.
    #[test]
    fn test_series_solve_power() {
        let mut circuit = create_series_circuit();

        circuit.compute_power();

        let edge_powers: Vec<f64> = circuit
            .graph
            .edge_weights()
            .map(|weight| weight.power.unwrap())
            .collect();
        let node_powers: Vec<f64> = circuit
            .graph
            .node_weights()
            .map(|weight| weight.power.unwrap())
            .collect();

        assert!(edge_powers[0] - 49.0 / 36.0 < 1e-10);
        assert!(edge_powers[1] - 25.0 / 36.0 < 1e-10);
        assert!(edge_powers[2] - 2.0 / 9.0 < 1e-10);
        assert!(edge_powers[3] - 2.0 / 9.0 < 1e-10);

        assert!(node_powers[0] - 7.0 / 3.0 < 1e-10);
        assert!(node_powers[1] - 5.0 / 4.0 < 1e-10);
        assert!(node_powers[2] - 4.0 / 3.0 < 1e-10);
        assert!(node_powers[3] - 0.0 < 1e-10);
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
        let source0 = VertexMetadata::new(Some(5.0), VertexTag::Source { tag: 0 });
        let source1 = VertexMetadata::new(Some(2.0), VertexTag::Source { tag: 1 });

        let v0 = VertexMetadata::new(None, VertexTag::Internal { tag: 2 });
        let v1 = VertexMetadata::new(None, VertexTag::Internal { tag: 3 });

        let sink = VertexMetadata::new(Some(0.0), VertexTag::Sink { tag: 4 });

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

    /// Test that the correct power values are found.
    #[test]
    fn test_multiple_solved_power() {
        let mut circuit = create_multiple_source_circuit();

        circuit.compute_power();

        let edge_powers: Vec<f64> = circuit
            .graph
            .edge_weights()
            .map(|weight| weight.power.unwrap())
            .collect();
        let node_powers: Vec<f64> = circuit
            .graph
            .node_weights()
            .map(|weight| weight.power.unwrap())
            .collect();

        assert!(edge_powers[0] - 676.0 / 121.0 < 1e-10);
        assert!(edge_powers[1] - 169.0 / 121.0 < 1e-10);
        assert!(edge_powers[2] - 169.0 / 121.0 < 1e-10);
        assert!(edge_powers[3] - 36.0 / 121.0 < 1e-10);
        assert!(edge_powers[4] - 256.0 / 121.0 < 1e-10);
        assert!(edge_powers[5] - 256.0 / 121.0 < 1e-10);

        assert!(node_powers[0] - 130.0 / 11.0 < 1e-10);
        assert!(node_powers[1] - 18.0 / 11.0 < 1e-10);
        assert!(node_powers[2] - 754.0 / 121.0 < 1e-10);
        assert!(node_powers[3] - 512.0 / 121.0 < 1e-10);
        assert!(node_powers[4] - 0.0 < 1e-10);
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
        let source = VertexMetadata::new(Some(c64::new(4.0, 0.0)), VertexTag::Source { tag: 0 });
        let v1 = VertexMetadata::new(None, VertexTag::Internal { tag: 1 });
        let v2 = VertexMetadata::new(None, VertexTag::Internal { tag: 2 });
        let sink = VertexMetadata::new(Some(c64::faer_zero()), VertexTag::Sink { tag: 3 });

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
            assert!((current - c64::new(8.0 / 5.0, 4.0 / 5.0)).abs() < 1e-10);
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

        assert!((solved_voltages[0] - c64::new(4.0, 0.0)).abs() < 1e-10);
        assert!((solved_voltages[1] - c64::new(12.0 / 5.0, 16.0 / 5.0)).abs() < 1e-10);
        assert!((solved_voltages[2] - c64::new(16.0 / 5.0, 8.0 / 5.0)).abs() < 1e-10);
        assert!((solved_voltages[3] - c64::faer_zero()).abs() < 1e-10);
    }

    /// Test that the correct power values are found.
    #[test]
    fn test_ac_power() {
        let mut circuit = create_ac_circuit();

        circuit.compute_power();

        let edge_powers: Vec<c64> = circuit
            .graph
            .edge_weights()
            .map(|weight| weight.power.unwrap())
            .collect();
        let node_powers: Vec<c64> = circuit
            .graph
            .node_weights()
            .map(|weight| weight.power.unwrap())
            .collect();

        assert!((edge_powers[0] - c64::new(128.0 / 25.0, -96.0 / 25.0)).abs() < 1e-10);
        assert!((edge_powers[1] - c64::new(-64.0 / 25.0, 48.0 / 25.0)).abs() < 1e-10);
        assert!((edge_powers[2] - c64::new(96.0 / 25.0, 128.0 / 25.0)).abs() < 1e-10);

        assert!((node_powers[0] - c64::new(32.0 / 5.0, 16.0 / 5.0)).abs() < 1e-10);
        assert!((node_powers[1] - c64::new(32.0 / 25.0, 176.0 / 25.0)).abs() < 1e-10);
        assert!((node_powers[2] - c64::new(96.0 / 25.0, 128.0 / 25.0)).abs() < 1e-10);
        assert!((node_powers[3] - c64::faer_zero()).abs() < 1e-10);
    }
}

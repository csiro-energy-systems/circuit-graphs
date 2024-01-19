use circuit_graphs::circuit_graph::*;
use criterion::{criterion_group, criterion_main, Criterion};
use faer_core::c64;
use faer_core::ComplexField;

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

    let e1 = EdgeMetadata::new(0, 1, EdgeType::new_component(c64::new(0.0, 0.5)));
    let e2 = EdgeMetadata::new(1, 2, EdgeType::new_component(c64::new(0.0, -1.0)));
    let e3 = EdgeMetadata::new(2, 3, EdgeType::new_component(c64::new(0.5, 0.0)));

    Circuit::new(vec![source, v1, v2, sink], vec![e1, e2, e3])
}

pub fn solve_ac(c: &mut Criterion) {
    c.bench_function("solve ac", |b| {
        b.iter_batched(
            || create_ac_circuit(),
            |mut circuit| circuit.solve_currents_and_voltages(),
            criterion::BatchSize::SmallInput,
        )
    });
}

/// Create a circuit with a transformer inside.
///
/// ```raw
///     _____
/// ----__R__---- ------|C(-----
/// |   2ohm    | |   -j2ohm   |
/// |+          $ $            $
/// V 6       20$T$30     j3ohmI
/// |-          $ $            $
/// |   _____   | |    _____   |
/// ----__R__---- -----__R__----
///     1ohm           2ohm
/// ```
fn create_transformer_circuit() -> Circuit<c64> {
    let source = VertexMetadata::new(Some(c64::new(6.0, 0.0)), 0, VertexType::Source);

    let transformer = VertexMetadata::new(None, 1, VertexType::TransformerSecondary);

    let internal1 = VertexMetadata::new(None, 2, VertexType::Internal);
    let internal2 = VertexMetadata::new(None, 3, VertexType::Internal);
    let internal3 = VertexMetadata::new(None, 4, VertexType::Internal);
    let internal4 = VertexMetadata::new(None, 5, VertexType::Internal);

    let primary_sink = VertexMetadata::new(Some(c64::faer_zero()), 6, VertexType::Sink);
    let secondary_sink = VertexMetadata::new(Some(c64::faer_zero()), 7, VertexType::Sink);

    let e0 = EdgeMetadata::new(0, 2, EdgeType::new_component(c64::new(0.5, 0.0)));
    let e1 = EdgeMetadata::new(2, 3, EdgeType::new_transformer(1, 20, 30));
    let e2 = EdgeMetadata::new(3, 6, EdgeType::new_component(c64::new(1.0, 0.0)));

    let e3 = EdgeMetadata::new(1, 4, EdgeType::new_component(c64::new(0.0, 0.5)));
    let e4 = EdgeMetadata::new(4, 5, EdgeType::new_component(c64::new(0.0, -1.0 / 3.0)));
    let e5 = EdgeMetadata::new(5, 7, EdgeType::new_component(c64::new(0.5, 0.0)));

    Circuit::new(
        vec![
            source,
            transformer,
            internal1,
            internal2,
            internal3,
            internal4,
            primary_sink,
            secondary_sink,
        ],
        vec![e0, e1, e2, e3, e4, e5],
    )
}

pub fn solve_transformer(c: &mut Criterion) {
    c.bench_function("solve transformer", |b| {
        b.iter_batched(
            || create_transformer_circuit(),
            |mut circuit| circuit.solve_currents_and_voltages(),
            criterion::BatchSize::SmallInput,
        )
    });
}

/// Set up a circuit with series transformers.
///
/// ```raw
/// ------------ ----------
/// |          | |        |
/// |          $ $       | |
/// |        1 $T$ 2     |R| 2ohm
/// |+         $ $       | |
/// V 17V      | |        |
/// |-         | ----------
/// |          | ----------
/// |          | |        |
/// |          $ $       | |
/// |        2 $T$ 3     |R| 1ohm
/// |          $ $       | |
/// |          | |        |
/// ------------ ----------
/// ```
fn create_series_transformer_circuit() -> Circuit<f64> {
    let source = VertexMetadata::new(Some(17.0), 0, VertexType::Source);
    let transformer0 = VertexMetadata::new(None, 1, VertexType::TransformerSecondary);
    let transformer1 = VertexMetadata::new(None, 2, VertexType::TransformerSecondary);
    let internal = VertexMetadata::new(None, 3, VertexType::Internal);
    let sink0 = VertexMetadata::new(Some(0.0), 4, VertexType::Sink);
    let sink1 = VertexMetadata::new(Some(0.0), 5, VertexType::Sink);
    let sink2 = VertexMetadata::new(Some(0.0), 6, VertexType::Sink);

    let e1 = EdgeMetadata::new(0, 3, EdgeType::new_transformer(1, 1, 2));
    let e2 = EdgeMetadata::new(3, 4, EdgeType::new_transformer(2, 2, 3));
    let e3 = EdgeMetadata::new(1, 5, EdgeType::new_component(0.5));
    let e4 = EdgeMetadata::new(2, 6, EdgeType::new_component(1.0));

    Circuit::new(
        vec![
            source,
            transformer0,
            transformer1,
            internal,
            sink0,
            sink1,
            sink2,
        ],
        vec![e1, e2, e3, e4],
    )
}

pub fn solve_series_transformer(c: &mut Criterion) {
    c.bench_function("solve series transformer", |b| {
        b.iter_batched(
            || create_series_transformer_circuit(),
            |mut circuit| circuit.solve_currents_and_voltages(),
            criterion::BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, solve_ac, solve_transformer, solve_series_transformer);
criterion_main!(benches);

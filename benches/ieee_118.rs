use std::fs::File;

use circuit_graphs::circuit_graph::*;
use criterion::{criterion_group, criterion_main, Criterion};
use faer_core::c64;
use faer_core::ComplexField;
use polars::prelude::*;

fn get_ieee_118_circuit() -> Circuit<c64> {
    let _unit_current: c64 = c64::new(9090.90909, 0.0);
    let unit_current: c64 = c64::new(11_000.0, 0.0);
    let _unit_power: c64 = c64::new(100_000.0, 0.0);

    let edge_file = File::open("ieee_118/data/transmission_line_data.parquet").unwrap();
    let edge_data = ParquetReader::new(edge_file)
        .finish()
        .unwrap()
        .select(["From Bus", "To Bus", "R (pu)", "X (pu)"])
        .unwrap();

    let mut edge_iters = edge_data
        .iter()
        .map(|series| series.iter())
        .collect::<Vec<_>>();

    let edges = (0..edge_data.height()).map(move |_| {
        EdgeMetadata::new(
            edge_iters[0].next().unwrap().try_extract().unwrap(),
            edge_iters[1].next().unwrap().try_extract().unwrap(),
            EdgeType::new_component(
                c64::new(
                    edge_iters[2].next().unwrap().try_extract().unwrap(),
                    edge_iters[3].next().unwrap().try_extract().unwrap(),
                )
                .faer_inv(),
            ),
        )
    });

    let generators_file = File::open("ieee_118/data/generator_data.parquet").unwrap();
    let generators_data = ParquetReader::new(generators_file)
        .finish()
        .unwrap()
        .select(["Bus No.", "Pmax (MW)", "Qmax (MVAR)"])
        .unwrap();

    let mut gen_iters = generators_data
        .iter()
        .map(|series| series.iter())
        .collect::<Vec<_>>();

    let mut source_nos = vec![];

    let sources = (0..generators_data.height())
        .map(|_| {
            let index = gen_iters[0].next().unwrap().try_extract().unwrap();
            source_nos.push(index);

            VertexMetadata::new(
                Some(
                    c64::new(
                        gen_iters[1].next().unwrap().try_extract().unwrap(),
                        gen_iters[2].next().unwrap().try_extract().unwrap(),
                    )
                    .faer_mul(unit_current.faer_inv()),
                ),
                index,
                VertexType::Source,
            )
        })
        .collect::<Vec<_>>();

    let nodes_file = File::open("ieee_118/data/bus_data.parquet").unwrap();
    let nodes_data = ParquetReader::new(nodes_file)
        .finish()
        .unwrap()
        .select(["Bus No."])
        .unwrap();

    let nodes = nodes_data
        .column("Bus No.")
        .unwrap()
        .iter()
        .map(|num| num.try_extract().unwrap())
        .filter(|index| !source_nos.contains(index))
        .map(|index| VertexMetadata::new(None, index, VertexType::Internal));

    let loads_file = File::open("ieee_118/data/bus_load_distribution_profile.parquet").unwrap();
    let loads_data = ParquetReader::new(loads_file).finish().unwrap();

    let mut load_iters = loads_data
        .iter()
        .map(|series| series.iter())
        .collect::<Vec<_>>();

    let loads = (0..loads_data.height()).map(move |_| {
        EdgeMetadata::new(
            load_iters[0].next().unwrap().try_extract().unwrap(),
            0,
            EdgeType::new_component(
                unit_current.faer_mul(unit_current).faer_mul(
                    c64::new(
                        load_iters[1].next().unwrap().try_extract().unwrap(),
                        load_iters[2].next().unwrap().try_extract().unwrap(),
                    )
                    .faer_inv(),
                ),
            ),
        )
    });

    Circuit::new(
        nodes
            .chain(vec![VertexMetadata::new(
                Some(c64::faer_zero()),
                0,
                VertexType::Sink,
            )])
            .chain(sources),
        edges.chain(loads),
    )
}

pub fn solve_ieee_118(c: &mut Criterion) {
    println!("Hey");
    c.bench_function("solve ieee", |b| {
        b.iter_batched(
            || get_ieee_118_circuit(),
            |mut circuit| circuit.solve_currents_and_voltages(),
            criterion::BatchSize::SmallInput,
        )
    });
}

criterion_group!(ieee, solve_ieee_118);
criterion_main!(ieee);

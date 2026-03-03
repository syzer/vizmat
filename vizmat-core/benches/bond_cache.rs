use std::path::Path;

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};

use vizmat_core::formats::parse_structure_by_extension;
use vizmat_core::structure::{infer_bonds_grid, resolve_bonds, BondInferenceSettings, Crystal};

fn load_structure_6vxx() -> Crystal {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../vizmat-app/assets/structures/proteins/6VXX.pdb");
    let contents = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read {path:?}: {e}"));
    parse_structure_by_extension("pdb", &contents)
        .unwrap_or_else(|e| panic!("failed to parse {path:?}: {e}"))
}

fn bond_inference_bench(c: &mut Criterion) {
    let crystal = load_structure_6vxx();
    let settings = BondInferenceSettings::default();

    let mut group = c.benchmark_group("bond_inference_6vxx");
    group.throughput(Throughput::Elements(crystal.atoms.len() as u64));

    group.bench_function("resolve_bonds_recompute", |b| {
        b.iter(|| {
            let (bonds, source) =
                resolve_bonds(black_box(&crystal), black_box(&settings));
            black_box((bonds.len(), source))
        });
    });

    group.finish();
}

fn cached_bond_bench(c: &mut Criterion) {
    let crystal = load_structure_6vxx();
    let settings = BondInferenceSettings::default();

    let (bonds, _source) = resolve_bonds(&crystal, &settings);

    let mut group = c.benchmark_group("bond_cache_6vxx");
    group.throughput(Throughput::Elements(crystal.atoms.len() as u64));

    group.bench_function("scan_cached_bonds", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for bond in black_box(&bonds) {
                total = total.saturating_add(bond.a).saturating_add(bond.b);
            }
            black_box(total)
        });
    });

    group.bench_function("infer_bonds_grid_direct", |b| {
        b.iter(|| {
            let bonds = infer_bonds_grid(black_box(&crystal), settings.tolerance_scale);
            black_box(bonds.len())
        });
    });

    group.finish();
}

criterion_group!(benches, bond_inference_bench, cached_bond_bench);
criterion_main!(benches);

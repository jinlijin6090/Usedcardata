import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ftc_enum import Network, make_species, species_from_smiles  # noqa: E402
from ftc_enum.rules import bimolecular, unimolecular  # noqa: E402


def S(smi):
    return species_from_smiles(smi)


def keyset(smiles_list):
    return frozenset(S(x).key for x in smiles_list)


def products(gen):
    out = set()
    for rtype, prods in gen:
        sps = [make_species(*p) for p in prods]
        if all(sps):
            out.add((rtype, frozenset(s.key for s in sps)))
    return out


def test_connectivity_key_distinguishes_h_count():
    assert S("C").key != S("[C]").key
    assert S("[CH3]").key != S("[CH2]").key


def test_bond_orders_and_free_valence():
    assert S("[C-]#[O+]").free_valence == 2       # CO 写成 [C]=O
    assert S("[CH3]").free_valence == 1
    assert S("[CH2]").free_valence == 2
    assert S("C=C").free_valence == 0 and S("C=C").n_pi == 1
    assert S("C#C").n_pi == 2
    assert S("CC[CH2]").smiles() == "[CH2]CC"
    assert S("[H][H]").nH == 2 and S("[H]").free_valence == 1


@pytest.mark.parametrize("a,b,prod", [
    ("[CH3]", "[CH3]", "CC"),            # 自由基结合
    ("[H]", "[C-]#[O+]", "[CH]=O"),      # H 加成到 CO
    ("[CH3]", "[C-]#[O+]", "C[C]=O"),    # CO 插入（乙酰基）
    ("[CH3]", "C=C", "CC[CH2]"),         # 自由基加成到 C=C
    ("[CH2]", "[CH2]", "C=C"),           # 卡宾二聚
    ("[H]", "[H]", "[H][H]"),
])
def test_b0f1(a, b, prod):
    sa, sb = S(a), S(b)
    got = products(bimolecular(sa, sb, b1f1=False))
    assert ("b0f1", keyset([prod])) in got


def test_b1f0_propyl_matches_yarp_sequential():
    # YARP sequential b1f0 对 CC[CH2] 给出的 5 组产物
    expected = {keyset(x) for x in (["C=C", "[CH3]"], ["C=CC", "[H]"], ["[CH2]", "[CH2]C"],
                                     ["[CH2]C[CH2]", "[H]"], ["[CH]CC", "[H]"])}
    got = {k for t, k in products(unimolecular(S("CC[CH2]"), b1f1=False))}
    assert got == expected


def test_b1f1_h_abstraction_and_disproportionation():
    got = products(bimolecular(S("[H][H]"), S("[CH3]"), b0f1=False))
    assert ("b1f1", keyset(["C", "[H]"])) in got
    got = products(bimolecular(S("C[CH2]"), S("C[CH2]"), b0f1=False))
    assert ("b1f1", keyset(["CC", "C=C"])) in got


def test_b1f1_closed_shell_acceptor_is_skipped_by_default():
    got = products(bimolecular(S("C"), S("C=C"), b0f1=False))
    assert not got
    got = products(bimolecular(S("C"), S("C=C"), b0f1=False, open_shell_acceptor=False))
    assert ("b1f1", keyset(["[CH3]", "C[CH2]"])) in got


def test_b1f1_unimolecular_keto_enol():
    got = products(unimolecular(S("CC=O"), b1f0=False))
    assert ("b1f1", keyset(["C=CO"])) in got


def test_network_two_heavy_atoms(tmp_path):
    import yaml
    with open(os.path.join(os.path.dirname(__file__), "..", "config.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    cfg["max_heavy"] = 2
    net = Network(cfg, str(tmp_path))
    for smi in cfg["seeds"]:
        net.add_seed(smi)
    summary = net.run(log_every=0)
    keys = set(net.species)
    for smi in ["C", "O", "CO", "C=O", "[CH]=O", "CC", "C=C", "C#C", "[CH2]C", "[OH]", "[H]"]:
        assert S(smi).key in keys, smi
    assert summary["n_species"] == 26
    assert summary["species_by_heavy_atoms"] == {0: 2, 1: 8, 2: 16}


def test_pack_unpack_preserves_atom_order():
    from ftc_enum.species import pack, unpack, species_from_key
    for smi in ["[H]", "[H][H]", "CCCCCC=C", "OCC(C)CC", "CCCCCCCCCC=CCCCCCCCCCCCO"]:
        s = S(smi)
        elems, hs, edges = unpack(pack(s))
        assert (tuple(elems), tuple(hs), tuple(edges)) == (s.elems, s.hs, s.edges)
        assert make_species(elems, hs, edges).key == s.key == species_from_key(s.key).key


def test_slimmed_partner_rehydrates_in_original_order(tmp_path):
    # 大分子先被处理并瘦身，之后的小分子以它为伙伴时，必须按原原子顺序重建
    import yaml
    with open(os.path.join(os.path.dirname(__file__), "..", "config.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    cfg["max_heavy"] = 7
    net = Network(cfg, str(tmp_path))
    seed = net.add_seed("CCCCCC=C")
    order = (seed.elems, seed.hs, seed.edges)
    net._mark_processed(seed)
    assert seed.nbrs is None
    full = net._rehydrate(seed)
    assert (full.elems, full.hs, full.edges) == order and full.id == seed.id and full.key == seed.key

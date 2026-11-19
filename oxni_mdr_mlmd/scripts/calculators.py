"""One entry point for the universal ML potentials compared in PLAN.md.

spec strings:
  mace:small | mace:medium | mace:medium-mpa-0 | mace:medium-omat-0 | mace:<path.model>
  mace:mace-matpes-pbe-0    (PBE without +U: the consistent starting point for Ni0 next to oxide)
  mace:mh-1@omat_pbe        (multi-head model trained on OMat/OMol/OC20/MatPES; "@head" picks the head)
  uma:uma-s-1p2p1/oc20      (fairchem-core >= 2; task: oc20 | oc22 | oc25 | omat | omol ...; needs HF access)
  sevennet:7net-omni/matpes_pbe (SevenNet; name/modal; older: 7net-mf-ompa/mpa)
  mattersim:MatterSim-v1.0.0-5M.pth
  orb:orb_v3_conservative_inf_omat
  dp:<path.pt>[@head]       (DeePMD-kit / DPA-3 fine-tuned model)
  emt                       (ASE EMT, pipeline tests only: no Mg/Si)

Only "mace:small" and "mace:medium" were run in this repository. The other MACE names are
those listed by the installed mace-torch 0.3.17; the remaining branches follow each
package's documented ASE interface - check them against the version you install.
"""
import os


def get_calculator(spec, device="cpu", dtype="float32"):
    kind, _, arg = spec.partition(":")
    if kind == "mace":
        if arg and os.path.exists(arg):
            from mace.calculators import MACECalculator
            return MACECalculator(model_paths=arg, device=device, default_dtype=dtype)
        from mace.calculators import mace_mp
        name, _, head = arg.partition("@")
        kw = {"head": head} if head else {}
        return mace_mp(model=name or "medium", device=device, default_dtype=dtype, **kw)
    if kind == "uma":
        from fairchem.core import FAIRChemCalculator, pretrained_mlip
        name, _, task = arg.partition("/")
        pu = pretrained_mlip.get_predict_unit(name or "uma-s-1p2p1", device=device)
        return FAIRChemCalculator(pu, task_name=task or "oc20")
    if kind == "sevennet":
        from sevenn.calculator import SevenNetCalculator
        name, _, modal = arg.partition("/")
        kw = {"modal": modal} if modal else {}
        return SevenNetCalculator(model=name or "7net-mf-ompa", device=device, **kw)
    if kind == "mattersim":
        from mattersim.forcefield import MatterSimCalculator
        return MatterSimCalculator(load_path=arg or "MatterSim-v1.0.0-5M.pth", device=device)
    if kind == "orb":
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.calculator import ORBCalculator
        model = getattr(pretrained, arg or "orb_v3_conservative_inf_omat")(device=device)
        return ORBCalculator(model, device=device)
    if kind == "dp":
        from deepmd.calculator import DP
        path, _, head = arg.partition("@")
        return DP(model=path, head=head) if head else DP(model=path)
    if kind == "emt":
        from ase.calculators.emt import EMT
        return EMT()
    raise ValueError("unknown calculator spec %r" % spec)

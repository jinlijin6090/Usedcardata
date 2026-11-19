"""One entry point for the universal ML potentials compared in PLAN.md.

spec strings:
  mace:small | mace:medium | mace:medium-mpa-0 | mace:medium-omat-0 | mace:<path.model>
  uma:uma-s-1p1/oc20        (fairchem-core >= 2; task: oc20 | omat | omol ...; needs HF access)
  sevennet:7net-mf-ompa/mpa (SevenNet; name/modal)
  mattersim:MatterSim-v1.0.0-5M.pth
  orb:orb_v3_conservative_inf_omat
  dp:<path.pt>[@head]       (DeePMD-kit / DPA-3 fine-tuned model)
  emt                       (ASE EMT, pipeline tests only: no Mg/Si)

Only the MACE and EMT branches were run in this repository; the others follow each
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
        return mace_mp(model=arg or "medium", device=device, default_dtype=dtype)
    if kind == "uma":
        from fairchem.core import FAIRChemCalculator, pretrained_mlip
        name, _, task = arg.partition("/")
        pu = pretrained_mlip.get_predict_unit(name or "uma-s-1p1", device=device)
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

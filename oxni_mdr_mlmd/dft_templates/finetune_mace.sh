#!/bin/bash
# Fine-tuning of a MACE foundation model on the OxNi DFT labels. Never run in this repository:
# check the flags against your mace version (mace_run_train --help). Train 3-4 seeds -> committee.
# E0s: isolated-atom energies computed with the SAME DFT settings (spin-polarised atoms in a box).
#
# FOUNDATION: default is mace-matpes-pbe-0 (PBE without +U), because the labels here are PBE-D3
# without U and the +U-trained models (small/medium, medium-mpa-0, medium-omat-0) put O/Ni on a
# different energy scale (see PLAN.md 2.2). Decide between candidates with benchmark_umlip.py.
# PT_TRAIN: mace-torch 0.3.17 does multi-head (replay) fine-tuning only for small/medium/large
# unless a file with pre-training configurations is given (--pt_train_file); without it the run
# falls back to plain fine-tuning, which forgets the foundation data faster. For
# mace-matpes-pbe-0 pass a subset of MatPES-PBE (extxyz) as PT_TRAIN.
set -e
FOUNDATION=${FOUNDATION:-mace-matpes-pbe-0}
PT_ARGS=""
if [ -n "${PT_TRAIN}" ]; then PT_ARGS="--multiheads_finetuning=True --pt_train_file=${PT_TRAIN}"; else PT_ARGS="--multiheads_finetuning=False"; fi
for SEED in 1 2 3 4; do
mace_run_train \
  --name="oxni_ft_s${SEED}" \
  --foundation_model="${FOUNDATION}" \
  ${PT_ARGS} \
  --train_file="train.extxyz" \
  --valid_file="valid.extxyz" \
  --test_file="test.extxyz" \
  --energy_key="REF_energy" \
  --forces_key="REF_forces" \
  --E0s="e0s.json" \
  --energy_weight=1.0 --forces_weight=10.0 \
  --lr=1e-4 --batch_size=4 --valid_batch_size=8 \
  --max_num_epochs=200 --patience=30 \
  --ema --ema_decay=0.99 --amsgrad \
  --default_dtype="float64" \
  --device=cuda \
  --seed=${SEED}
done

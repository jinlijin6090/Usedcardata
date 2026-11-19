#!/bin/bash
# Multi-head fine-tuning of a MACE foundation model on the OxNi DFT labels.
# Check flags against your mace version (mace_run_train --help). Train 3-4 seeds -> committee.
# E0s: isolated-atom energies computed with the SAME DFT settings (spin-polarised atoms in a box).
set -e
FOUNDATION=${FOUNDATION:-medium-omat-0}   # or medium-mpa-0 / small; decided by benchmark_umlip.py
for SEED in 1 2 3 4; do
mace_run_train \
  --name="oxni_ft_s${SEED}" \
  --foundation_model="${FOUNDATION}" \
  --multiheads_finetuning=True \
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

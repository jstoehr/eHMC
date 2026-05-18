#!/bin/bash

# --------------------------------------------------
# --- Default settings
# --------------------------------------------------

MODELS="MVNorm banana BLP"
SEEDS="1 2 3 4 5"

NTRAIN="50000"
ESS_TARGETS="2500 5000"
N_WARMUPS="2000 5000"
RESAMPLING=1
MAX_TRIES=100
TRANSFORMS=5
HIDDEN_FEATURES="(32,32)"
BATCH_SIZE=2048
LEARNING_RATE=1e-3

REPO_ROOT="."
DEVICE="cpu"

# --------------------------------------------------
# --- Parse arguments
# --------------------------------------------------

while [[ $# -gt 0 ]]; do
    case $1 in
        --model)
            MODELS="$2"
            shift 2
            ;;
        --seed)
            SEEDS="$2"
            shift 2
            ;;
        --n_train)
            NTRAIN="$2"
            shift 2
            ;;
        --ess_target)
            ESS_TARGETS="$2"
            shift 2
            ;;
        --n_warmup)
            N_WARMUPS="$2"
            shift 2
            ;;
        --resampling)
            RESAMPLING="$2"
            shift 2
            ;;
        --max_tries)
            MAX_TRIES="$2"
            shift 2
            ;;
        --transforms)
            TRANSFORMS="$2"
            shift 2
            ;;
        --hidden_features)
            HIDDEN_FEATURES="$2"
            shift 2
            ;;
         --batch_size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --learning_rate)
            LEARNING_RATE="$2"
            shift 2
            ;;
        --repo_root)
            REPO_ROOT="$2"
            shift 2
            ;;
        --device)
            DEVICE="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# --------------------------------------------------
# --- Run experiments
# --------------------------------------------------

if [[ "$SEEDS" == *:* ]]; then
    IFS=':' read -r SEED_START SEED_END <<< "$SEEDS"
    SEED_LIST=$(seq "$SEED_START" "$SEED_END")
else
    SEED_LIST="$SEEDS"
fi

for model in $MODELS; do
    for seed in $SEED_LIST; do
        for ess_target in $ESS_TARGETS; do
            for n_warmup in $N_WARMUPS; do

                echo "--------------------------------------------------------"
                echo "Model: $model | Seed: $seed | n_train: $NTRAIN | ESS_TARGET: $ess_target | n_warmup: $n_warmup"
                echo "--------------------------------------------------------"

                python examples/scripts/run_calibration.py \
                    --repo_root "$REPO_ROOT" \
                    --model "$model" \
                    --seed "$seed" \
                    --n_train "$NTRAIN" \
                    --ess_target "$ess_target" \
                    --n_warmup "$n_warmup" \
                    --resampling "$RESAMPLING" \
                    --max_tries "$MAX_TRIES" \
                    --transforms "$TRANSFORMS" \
                    --hidden_features "$HIDDEN_FEATURES" \
                    --batch_size "$BATCH_SIZE" \
                    --learning_rate "$LEARNING_RATE" \
                    --device "$DEVICE"
            done    
        done
    done
done


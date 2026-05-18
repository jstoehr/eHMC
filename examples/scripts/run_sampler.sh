#!/bin/bash

# --------------------------------------------------
# --- Default settings
# --------------------------------------------------

MODELS="MVNorm banana BLP"
ALGOS="ehmc prhmc nuts"
SEEDS="1:20"

CHAINS=50
WARMUPS="2000 5000"
ITER=2000
DELTAS="0.651 0.8"
METRIC="diag"

NTRAIN=50000
ESS_TARGETS="2500 5000"
ESS_TRAIN=0.8
RESAMPLING=1
EPOCH=1
REFRESH=0.75
RAND=0

REPO_ROOT="."

# --------------------------------------------------
# --- Parse arguments
# --------------------------------------------------

while [[ $# -gt 0 ]]; do
    case $1 in
        --model)
            MODELS="$2"
            shift 2
            ;;
        --algo)
            ALGOS="$2"
            shift 2
            ;;
        --seed)
            SEEDS="$2"
            shift 2
            ;;
        --chains)
            CHAINS="$2"
            shift 2
            ;;
        --warmup)
            WARMUPS="$2"
            shift 2
            ;;
        --iter)
            ITER="$2"
            shift 2
            ;;
        --delta)
            DELTAS="$2"
            shift 2
            ;;
        --metric)
            METRIC="$2"
            shift 2
            ;;
        --n_train)
            NTRAIN="$2"
            shift 2
            ;;
        --ess_train)
            ESS_TRAIN="$2"
            shift 2
            ;;
        --ess_target)
            ESS_TARGETS="$2"
            shift 2
            ;;
        --resampling)
            RESAMPLING="$2"
            shift 2
            ;;
        --epoch)
            EPOCH="$2"
            shift 2
            ;;
        --refresh)
            REFRESH="$2"
            shift 2
            ;;
        --rand)
            RAND="$2"
            shift 2
            ;;
        --repo_root)
            REPO_ROOT="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# --------------------------------------------------
# --- Expand seed range
# --------------------------------------------------

if [[ "$SEEDS" == *:* ]]; then
    IFS=':' read -r SEED_START SEED_END <<< "$SEEDS"
    SEED_LIST=$(seq "$SEED_START" "$SEED_END")
else
    SEED_LIST="$SEEDS"
fi

# --------------------------------------------------
# --- Run experiments
# --------------------------------------------------

for model in $MODELS; do
    for algo in $ALGOS; do
        for seed in $SEED_LIST; do
            for warmup in $WARMUPS; do
                for delta in $DELTAS; do

                    if [[ "$algo" == "nuts" ]]; then
                        echo "--------------------------------------------------------"
                        echo "Algorithm: $algo | Model: $model | Seed: $seed | warmup: $warmup | delta: $delta"
                        echo "--------------------------------------------------------"

                        Rscript examples/scripts/run_sampler.R \
                            --repo_root "$REPO_ROOT" \
                            --model "$model" \
                            --algo "$algo" \
                            --seed "$seed" \
                            --chains "$CHAINS" \
                            --warmup "$warmup" \
                            --iter "$ITER" \
                            --delta "$delta" \
                            --metric "$METRIC"

                    else
                        for ess_target in $ESS_TARGETS; do
                            echo "--------------------------------------------------------"
                            echo "Algorithm: $algo | Model: $model | Seed: $seed | ESS: $ess_target | warmup: $warmup | delta: $delta"
                            echo "--------------------------------------------------------"

                            Rscript examples/scripts/run_sampler.R \
                                --repo_root "$REPO_ROOT" \
                                --model "$model" \
                                --algo "$algo" \
                                --seed "$seed" \
                                --chains "$CHAINS" \
                                --warmup "$warmup" \
                                --iter "$ITER" \
                                --delta "$delta" \
                                --metric "$METRIC" \
                                --n_train "$NTRAIN" \
                                --ess_train "$ESS_TRAIN" \
                                --ess_target "$ess_target" \
                                --resampling "$RESAMPLING" \
                                --epoch "$EPOCH" \
                                --refresh "$REFRESH" \
                                --rand "$RAND"
                        done
                    fi

                done
            done
        done
    done
done
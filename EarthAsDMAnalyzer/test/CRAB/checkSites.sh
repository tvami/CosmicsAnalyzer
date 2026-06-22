#!/bin/bash

datasets=(
"/Cosmics/Run2022B-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2022C-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2022D-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2022D-CosmicTP-PromptReco-v2/RAW-RECO"
"/Cosmics/Run2022D-CosmicTP-PromptReco-v3/RAW-RECO"
"/Cosmics/Run2022E-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2022F-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2022G-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2023A-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2023A-CosmicTP-PromptReco-v2/RAW-RECO"
"/Cosmics/Run2023B-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2023C-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2023C-CosmicTP-PromptReco-v2/RAW-RECO"
"/Cosmics/Run2023C-CosmicTP-PromptReco-v3/RAW-RECO"
"/Cosmics/Run2023C-CosmicTP-PromptReco-v4/RAW-RECO"
"/Cosmics/Run2023D-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2023D-CosmicTP-PromptReco-v2/RAW-RECO"
"/Cosmics/Run2023E-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2023F-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2024A-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2024B-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2024C-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2024D-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2024E-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2024E-CosmicTP-PromptReco-v2/RAW-RECO"
"/Cosmics/Run2024F-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2024G-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2024H-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2024I-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2024I-CosmicTP-PromptReco-v2/RAW-RECO"
"/Cosmics/Run2024J-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2025A-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2025A-CosmicTP-PromptReco-v2/RAW-RECO"
"/Cosmics/Run2025B-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2025C-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2025C-CosmicTP-PromptReco-v2/RAW-RECO"
"/Cosmics/Run2025D-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2025E-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2025F-CosmicTP-PromptReco-v1/RAW-RECO"
"/Cosmics/Run2025F-CosmicTP-PromptReco-v2/RAW-RECO"
"/Cosmics/Run2025G-CosmicTP-PromptReco-v1/RAW-RECO"
)

for ds in "${datasets[@]}"; do
    echo "=============================="
    echo "Dataset: $ds"
    echo "Sites:"
    dasgoclient --query="site dataset=$ds" | sort -u
    echo
done


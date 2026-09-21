#!/usr/bin/env python3
"""
Multi-Agent Therapeutic lncRNA Pre-Screening Pipeline
------------------------------------------------------
Simulates an agentic workflow integrating lncRank-Spec mathematical engines
(Tau Index & lncScore), DepMap CRISPR dependency scoring, lncATLAS subcellular
localization, and cross-species orthology mapping for Glioblastoma (GBM) and
Triple-Negative Breast Cancer (TNBC).
"""

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class TranscriptData:
    gene_symbol: str
    ensembl_id: str
    disease_context: str
    target_expression_tpm: float
    healthy_organ_expressions: Dict[str, float]  # organ -> TPM
    depmap_dependency_score: float               # Chronos/CRISPR score (< -0.5 indicates dependency)
    lncatlas_rci: float                          # Nuclear vs Cytoplasmic RCI (> 0 = Nuclear, < 0 = Cytoplasmic)
    mouse_tau: Optional[float]                   # Mouse ortholog Tau index
    rat_tau: Optional[float]                     # Rat ortholog Tau index


@dataclass
class PrescreeningResult:
    gene_symbol: str
    disease_context: str
    tau_index: float
    lnc_score: float
    depmap_status: str
    subcellular_location: str
    recommended_modality: str
    safety_status: str
    species_conservation: str
    overall_rank_score: float


class IngestionFilterAgent:
    """Agent 1: Ingests expression data and applies dynamic detection thresholds."""
    
    def __init__(self, tpm_threshold: float = 1.5):
        self.tpm_threshold = tpm_threshold

    def filter_noise(self, dataset: List[TranscriptData]) -> List[TranscriptData]:
        filtered = []
        for item in dataset:
            if item.target_expression_tpm >= self.tpm_threshold:
                filtered.append(item)
            else:
                print(f"[IngestionAgent] Filtered out {item.gene_symbol}: Target TPM ({item.target_expression_tpm}) below threshold ({self.tpm_threshold})")
        return filtered


class RankingMathAgent:
    """Agent 2: Calculates Tau specificity index and lncScore."""
    
    @staticmethod
    def calculate_tau(target_tpm: float, healthy_exprs: Dict[str, float]) -> float:
        """
        Calculates Yanai's Tau index across all evaluated tissues.
        Tau = sum(1 - x_i / x_max) / (N - 1)
        """
        all_exprs = list(healthy_exprs.values()) + [target_tpm]
        x_max = max(all_exprs)
        if x_max == 0:
            return 0.0
        n = len(all_exprs)
        tau = sum(1.0 - (x / x_max) for x in all_exprs) / (n - 1)
        return round(min(max(tau, 0.0), 1.0), 4)

    @staticmethod
    def calculate_lnc_score(tau: float, mean_target_tpm: float) -> float:
        """
        lncScore = Tau^2 * log2(Mean Expression + 1)
        """
        return round((tau ** 2) * math.log2(mean_target_tpm + 1.0), 2)


class VulnerabilityAgent:
    """Agent 3: Evaluates DepMap contextual dependency & lncATLAS subcellular localization."""

    @staticmethod
    def evaluate(depmap_score: float, lncatlas_rci: float) -> tuple[str, str, str]:
        # DepMap interpretation
        if depmap_score <= -0.6:
            dep_status = "Critical Dependency"
        elif depmap_score <= -0.3:
            dep_status = "Moderate Dependency"
        else:
            dep_status = "Non-Essential"

        # lncATLAS Subcellular Localization & Therapeutic Modality
        if lncatlas_rci >= 1.0:
            location = "Nuclear Enriched"
            modality = "Antisense Oligonucleotide (ASO) / GapmeR"
        elif lncatlas_rci <= -1.0:
            location = "Cytoplasmic Enriched"
            modality = "siRNA / miRNA Sponge Disrupter"
        else:
            location = "Biphasic / Pervasive"
            modality = "Dual ASO / Small Molecule"

        return dep_status, location, modality


class SafetyOrthologyAgent:
    """Agent 4: Assesses off-target organ safety and cross-species conservation."""

    @staticmethod
    def assess_safety(healthy_exprs: Dict[str, float], threshold: float = 1.0) -> str:
        leaky_organs = [organ for organ, tpm in healthy_exprs.items() if tpm > threshold]
        if not leaky_organs:
            return "Zero Off-Target Risk (Strictly Restricted)"
        else:
            return f"Caution: Expression in {', '.join(leaky_organs)}"

    @staticmethod
    def assess_conservation(mouse_tau: Optional[float], rat_tau: Optional[float]) -> str:
        if mouse_tau is not None and rat_tau is not None:
            if mouse_tau >= 0.8 and rat_tau >= 0.8:
                return "High Conservation (Mouse & Rat Validated)"
            elif mouse_tau >= 0.7 or rat_tau >= 0.7:
                return "Moderate Conservation"
        elif mouse_tau is not None and mouse_tau >= 0.8:
            return "Mouse Conserved"
        return "Primate / Human Specific"


class MultiAgentPrescreeningOrchestrator:
    """Orchestrates all agents to run an automated pre-screening pipeline."""

    def __init__(self, tpm_noise_cutoff: float = 1.5):
        self.ingestion_agent = IngestionFilterAgent(tpm_threshold=tpm_noise_cutoff)
        self.math_agent = RankingMathAgent()
        self.vulnerability_agent = VulnerabilityAgent()
        self.safety_agent = SafetyOrthologyAgent()

    def run_pipeline(self, raw_candidates: List[TranscriptData]) -> List[PrescreeningResult]:
        print("=" * 80)
        print("STARTING MULTI-AGENT THERAPEUTIC PRE-SCREENING PIPELINE")
        print("=" * 80)
        
        # Step 1: Noise Filtering
        clean_candidates = self.ingestion_agent.filter_noise(raw_candidates)
        results = []

        # Process each candidate through agents
        for cand in clean_candidates:
            # Step 2: Specificity & Ranking Math
            tau = self.math_agent.calculate_tau(cand.target_expression_tpm, cand.healthy_organ_expressions)
            lnc_score = self.math_agent.calculate_lnc_score(tau, cand.target_expression_tpm)

            # Step 3: Vulnerability Evaluation
            dep_status, location, modality = self.vulnerability_agent.evaluate(
                cand.depmap_dependency_score, cand.lncatlas_rci
            )

            # Step 4: Safety & Orthology
            safety_status = self.safety_agent.assess_safety(cand.healthy_organ_expressions)
            conservation = self.safety_agent.assess_conservation(cand.mouse_tau, cand.rat_tau)

            # Overall Composite Rank Score (lncScore weighted by dependency & safety)
            dep_multiplier = 1.3 if "Critical" in dep_status else (1.1 if "Moderate" in dep_status else 0.8)
            safety_multiplier = 1.2 if "Zero" in safety_status else 0.7
            overall_score = round(lnc_score * dep_multiplier * safety_multiplier, 2)

            res = PrescreeningResult(
                gene_symbol=cand.gene_symbol,
                disease_context=cand.disease_context,
                tau_index=tau,
                lnc_score=lnc_score,
                depmap_status=dep_status,
                subcellular_location=location,
                recommended_modality=modality,
                safety_status=safety_status,
                species_conservation=conservation,
                overall_rank_score=overall_score
            )
            results.append(res)

        # Sort results by overall composite rank score
        results.sort(key=lambda x: x.overall_rank_score, reverse=True)
        return results


def print_report(results: List[PrescreeningResult]):
    """Prints a clean formatted CLI report of the screening findings."""
    print("\n" + "#" * 80)
    print("FINAL PRE-SCREENING TARGET PRIORITIZATION REPORT")
    print("#" * 80)
    
    contexts = set(r.disease_context for r in results)
    for ctx in contexts:
        print(f"\n>>> DISEASE CONTEXT: {ctx.upper()} <<<")
        print("-" * 80)
        ctx_results = [r for r in results if r.disease_context == ctx]
        for rank, r in enumerate(ctx_results, start=1):
            print(f"Rank #{rank}: {r.gene_symbol}")
            print(f"  * Tau Index (τ)        : {r.tau_index:.4f}")
            print(f"  * lncScore             : {r.lnc_score}")
            print(f"  * Overall Rank Score   : {r.overall_rank_score}")
            print(f"  * DepMap Vulnerability : {r.depmap_status}")
            print(f"  * Subcellular Location : {r.subcellular_location}")
            print(f"  * Recommended Modality : {r.recommended_modality}")
            print(f"  * Safety Profile       : {r.safety_status}")
            print(f"  * Cross-Species Status : {r.species_conservation}")
            print("-" * 80)


if __name__ == "__main__":
    # Sample Dataset for Glioblastoma (GBM) and Triple-Negative Breast Cancer (TNBC)
    sample_transcripts = [
        # Glioblastoma Candidates
        TranscriptData(
            gene_symbol="lncGRS-1",
            ensembl_id="ENSG00000282100",
            disease_context="Glioblastoma",
            target_expression_tpm=28.5,
            healthy_organ_expressions={"Heart": 0.05, "Liver": 0.1, "Kidney": 0.0, "Healthy_Brain": 0.2},
            depmap_dependency_score=-0.82,
            lncatlas_rci=2.1,  # Nuclear Enriched
            mouse_tau=None,     # Primate/Human specific
            rat_tau=None
        ),
        TranscriptData(
            gene_symbol="FOXM1-AS",
            ensembl_id="ENSG00000250110",
            disease_context="Glioblastoma",
            target_expression_tpm=18.2,
            healthy_organ_expressions={"Heart": 0.2, "Liver": 0.1, "Kidney": 0.1, "Healthy_Brain": 0.3},
            depmap_dependency_score=-0.68,
            lncatlas_rci=1.8,  # Nuclear
            mouse_tau=0.88,
            rat_tau=0.85
        ),
        TranscriptData(
            gene_symbol="HIF1A-AS2",
            ensembl_id="ENSG00000237290",
            disease_context="Glioblastoma",
            target_expression_tpm=14.1,
            healthy_organ_expressions={"Heart": 0.4, "Liver": 0.2, "Kidney": 0.3, "Healthy_Brain": 0.5},
            depmap_dependency_score=-0.54,
            lncatlas_rci=1.2,
            mouse_tau=0.82,
            rat_tau=0.80
        ),
        
        # Triple-Negative Breast Cancer Candidates
        TranscriptData(
            gene_symbol="TROJAN",
            ensembl_id="ENSG00000261100",
            disease_context="TNBC",
            target_expression_tpm=22.4,
            healthy_organ_expressions={"Heart": 0.0, "Liver": 0.05, "Kidney": 0.0, "Healthy_Breast": 0.1},
            depmap_dependency_score=-0.76,
            lncatlas_rci=1.9,  # Nuclear
            mouse_tau=None,
            rat_tau=None
        ),
        TranscriptData(
            gene_symbol="LASTR",
            ensembl_id="ENSG00000265700",
            disease_context="TNBC",
            target_expression_tpm=19.8,
            healthy_organ_expressions={"Heart": 0.1, "Liver": 0.1, "Kidney": 0.2, "Healthy_Breast": 0.2},
            depmap_dependency_score=-0.65,
            lncatlas_rci=1.4,
            mouse_tau=0.84,
            rat_tau=0.81
        ),
        TranscriptData(
            gene_symbol="GATA3-AS1",
            ensembl_id="ENSG00000242500",
            disease_context="TNBC",
            target_expression_tpm=16.5,
            healthy_organ_expressions={"Heart": 0.1, "Liver": 0.3, "Kidney": 0.1, "Healthy_Breast": 0.4},
            depmap_dependency_score=-0.58,
            lncatlas_rci=-1.5, # Cytoplasmic
            mouse_tau=0.79,
            rat_tau=0.77
        ),
        TranscriptData(
            gene_symbol="LINK-A",
            ensembl_id="ENSG00000211390",
            disease_context="TNBC",
            target_expression_tpm=25.1,
            healthy_organ_expressions={"Heart": 0.2, "Liver": 0.1, "Kidney": 0.1, "Healthy_Breast": 0.3},
            depmap_dependency_score=-0.71,
            lncatlas_rci=-1.8, # Cytoplasmic
            mouse_tau=0.86,
            rat_tau=0.83
        ),
        # Noise Control Example
        TranscriptData(
            gene_symbol="NOISE-RNA1",
            ensembl_id="ENSG00000999999",
            disease_context="Glioblastoma",
            target_expression_tpm=0.4, # Below 1.5 TPM cutoff
            healthy_organ_expressions={"Heart": 0.1, "Liver": 0.1, "Kidney": 0.1, "Healthy_Brain": 0.1},
            depmap_dependency_score=-0.1,
            lncatlas_rci=0.0,
            mouse_tau=0.1,
            rat_tau=0.1
        )
    ]

    orchestrator = MultiAgentPrescreeningOrchestrator(tpm_noise_cutoff=1.5)
    final_results = orchestrator.run_pipeline(sample_transcripts)
    print_report(final_results)

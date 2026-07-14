from __future__ import annotations

from collections import Counter
from dataclasses import replace

from .models import ClusterInfo, MappingDecision
from .mmseqs import ClusterIndex


def connect_proteins_to_clusters(
    decisions: list[MappingDecision], cluster_index: ClusterIndex
) -> list[MappingDecision]:
    cluster_by_member = cluster_index.clusters_for({
        decision.uniref90_id for decision in decisions
        if decision.status == "mapped" and decision.uniref90_id
    })
    connected: list[MappingDecision] = []
    for decision in decisions:
        cluster_id = ""
        status = decision.status
        detail = decision.detail
        if decision.status == "mapped" and decision.uniref90_id:
            cluster_id = cluster_by_member.get(decision.uniref90_id, "")
            if not cluster_id:
                status = "missing-mmseqs-assignment"
                detail = "mapped UniRef90 identifier has no MMseqs2 cluster assignment"
        connected.append(replace(
            decision, mmseqs_cluster_id=cluster_id, status=status, detail=detail,
        ))
    return connected


def retained_cluster_info(
    decisions: list[MappingDecision], cluster_index: ClusterIndex
) -> dict[str, ClusterInfo]:
    labelled_pairs = {
        (decision.mmseqs_cluster_id, decision.protein_id) for decision in decisions
        if decision.status == "mapped" and decision.mmseqs_cluster_id
    }
    labelled = Counter(cluster_id for cluster_id, _ in labelled_pairs)
    retained_ids = set(labelled)
    sizes = cluster_index.sizes_for(retained_ids)
    retained: dict[str, ClusterInfo] = {}
    for cluster_id in sorted(labelled):
        retained[cluster_id] = ClusterInfo(
            cluster_id=cluster_id,
            member_count=sizes[cluster_id],
            labelled_protein_count=labelled[cluster_id],
        )
    return retained


def members_for_retained(
    cluster_index: ClusterIndex, retained: dict[str, ClusterInfo]
) -> list[tuple[str, str]]:
    retained_ids = set(retained)
    return [
        (cluster_id, member_id)
        for cluster_id, member_id in cluster_index.iter_assignments()
        if cluster_id in retained_ids
    ]


def mapping_counters(decisions: list[MappingDecision]) -> dict[str, int]:
    counts = Counter(decision.status for decision in decisions)
    counts["canonical_sequence_available"] = sum(
        decision.canonical_sequence_available for decision in decisions
    )
    counts["mapped_to_cluster"] = sum(bool(decision.mmseqs_cluster_id) for decision in decisions)
    return {key: int(counts[key]) for key in sorted(counts)}

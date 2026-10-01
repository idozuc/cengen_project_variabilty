"""Shared neuron-family and ligand/receptor reference definitions."""

FAMILIES = {
    "motor_body": ["DA", "DB", "VA", "VB", "AS", "VD_DD", "VC", "VC_4_5"],
    "motor_head": ["SMB", "SMD", "SMD_stressed", "RMD_DV", "RMD_LR", "RME_DV", "RME_LR", "RMF", "RMG", "RMH"],
    "sensory": ["AFD", "AWA", "AWB", "AWC_ON", "AWC_OFF", "ASEL", "ASER", "ASH", "ASI", "ASJ", "ASK", "ADF", "ADL", "AQR", "BAG", "FLP", "PDE", "URX", "URY", "PHA", "PHB", "PHC", "OLQ", "CEP", "IL1", "IL2_DV", "IL2_LR"],
    "interneuron": ["RIA", "RIB", "RIC", "RID", "RIF", "RIG", "RIH", "RIM", "RIP", "RIR", "RIS", "RIV", "RIV_stressed", "AIB", "AIA", "AIM", "AIN", "AIY", "AVA", "AVB", "AVD", "AVE", "AVF", "AVG", "AVH", "AVJ", "AVK", "AVL", "DVA", "DVB", "DVC", "LUA", "PVC", "PVD", "PVM", "PVP", "PVQ", "PVR", "PVW", "PLM", "PLN", "BDU", "ADA", "ADE", "AUA", "SAA", "SAB", "SDQ", "URA", "URB", "AVM", "ALM", "PQR"],
    "pharyngeal": ["I1", "I4", "I5", "I6", "M1", "M2", "M3", "M5", "MI", "NSM", "MC"],
    "other": ["HSN", "SIA", "SIB"],
}

CELL_TYPE_TO_FAMILY = {cell_type: family for family, types in FAMILIES.items() for cell_type in types}

LIGAND_TO_RECEPTOR = {
    "flp-1": ["npr-22"], "flp-2": ["npr-4", "npr-11"], "flp-5": ["npr-22"],
    "flp-7": ["ckr-2"], "flp-8": ["frpr-6", "frpr-11"], "flp-9": ["frpr-8"],
    "flp-10": ["frpr-10"], "flp-11": ["frpr-3"], "flp-12": ["frpr-12"],
    "flp-13": ["frpr-4", "frpr-5"], "flp-14": ["frpr-13"],
    "flp-15": ["frpr-1", "frpr-11"], "flp-17": ["frpr-2"],
    "flp-18": ["npr-1", "npr-4", "npr-11"], "flp-19": ["npr-12"],
    "flp-20": ["frpr-4"], "flp-21": ["frpr-3"], "flp-22": ["frpr-4", "frpr-5"],
    "flp-32": ["frpr-12"], "nlp-1": ["npr-11"], "nlp-3": ["npr-22"],
    "nlp-12": ["ckr-2"], "nlp-14": ["ckr-1"], "nlp-15": ["npr-10", "npr-12"],
    "nlp-21": ["npr-1"], "nlp-29": ["frpr-4"], "pdf-1": ["pdfr-1"],
    "pdf-2": ["pdfr-1"],
}

RECEPTOR_TO_LIGANDS = {}
for ligand, receptors in LIGAND_TO_RECEPTOR.items():
    for receptor in receptors:
        RECEPTOR_TO_LIGANDS.setdefault(receptor, []).append(ligand)

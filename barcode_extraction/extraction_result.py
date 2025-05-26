class DetectedElement:
    def __init__(self, start=-1, end=-1, score=-1, seq=None):
        self.start = start
        self.end = end
        self.score = score
        self.seq = seq


class ExtractionResult:
    def __init__(self, read_id: str, strand: str, detected_results: dict):
        self.read_id = read_id
        self.stand = strand
        self.detected_results = detected_results

    # TODO: add simple serialziation

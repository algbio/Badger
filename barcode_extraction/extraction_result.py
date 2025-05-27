from collections import defaultdict


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


class ReadStats:
    def __init__(self):
        self.read_count = 0
        self.pattern_counts = defaultdict(int)

    def add_read(self, barcode_detection_result):
        self.read_count += 1
        for el in barcode_detection_result.detected_results:
            if barcode_detection_result.detected_results[el].start != -1:
                self.pattern_counts[el] += 1

    def __str__(self):
        human_readable_str =  ("Total reads:\t%d\n" % self.read_count)
        for a in self.pattern_counts:
            human_readable_str += "%s:\t%d\n" % (a, self.pattern_counts[a])
        return human_readable_str

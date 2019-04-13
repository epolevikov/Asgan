from asgan.output_generator import pretty_number


class PafHit:
    def __init__(self, query_name, query_len, query_start, query_end,
                 target_name, target_len, target_start, target_end):
        self.id = None

        self.query_name = query_name
        self.query_len = query_len
        self.query_start = query_start
        self.query_end = query_end

        self.target_name = target_name
        self.target_len = target_len
        self.target_start = target_start
        self.target_end = target_end

    def signed_id(self):
        return ["-", "+"][self.id > 0] + str(abs(self.id))

    def query_hit_length(self):
        return self.query_end - self.query_start

    def target_hit_length(self):
        return self.target_end - self.target_start

    def __str__(self):
        query_info = "{}\t{}\t{}\t{}\t".format(
            self.query_name,
            pretty_number(self.query_len),
            pretty_number(self.query_start),
            pretty_number(self.query_end))
        target_info = "{}\t{}\t{}\t{}".format(
            self.target_name,
            pretty_number(self.target_len),
            pretty_number(self.target_start),
            pretty_number(self.target_end))

        return query_info + target_info


def process_raw_hits(raw_hits):
    raw_hits = _filter_hits_by_len(raw_hits)

    processed_hits = []
    for raw_hit in raw_hits:
        processed_hit = _process_raw_hit(raw_hit)
        processed_hits.append(processed_hit)

    united_hits = _unite_processed_hits(processed_hits)

    for i, hit in enumerate(united_hits):
        hit.id = i + 1

    hits = []
    for hit in united_hits:
        hits.append(hit)
        hits.append(_complement_hit(hit))

    return hits


def filter_repeats(raw_hits, repeats_query, repeats_target):
    filtered_hits = []

    for hit in raw_hits:
        if hit.query_name in repeats_query:
            continue

        if hit.target_name in repeats_target:
            continue

        filtered_hits.append(hit)

    return filtered_hits


def _filter_hits_by_len(raw_hits, min_hit_length=50000):
    filtered_hits = []

    for raw_hit in raw_hits:
        if raw_hit.query_hit_length() < min_hit_length:
            continue

        if raw_hit.target_hit_length() < min_hit_length:
            continue

        filtered_hits.append(raw_hit)

    return filtered_hits


def _process_raw_hit(raw_hit):
    query_name = raw_hit.query_name + "+"
    query_start = raw_hit.query_start
    query_end = raw_hit.query_end
    target_name = raw_hit.target_name + raw_hit.strand

    if raw_hit.strand == "+":
        target_start = raw_hit.target_start
        target_end = raw_hit.target_end
    else:
        target_start = raw_hit.target_len - raw_hit.target_end
        target_end = raw_hit.target_len - raw_hit.target_start

    processed_hit = PafHit(
        query_name, raw_hit.query_len, query_start, query_end,
        target_name, raw_hit.target_len, target_start, target_end)

    return processed_hit


def _unite_processed_hits(processed_hits):
    processed_hits.sort(key=lambda hit: (-hit.query_len, hit.query_name,
                                         hit.query_start))

    united_hits = []
    curr_hit = processed_hits[0]
    max_hits_dist = 5 * 10**4

    for i in range(1, len(processed_hits)):
        next_hit = processed_hits[i]

        query_hits_dist = next_hit.query_start - curr_hit.query_end
        target_hits_dist = next_hit.target_start - curr_hit.target_end

        if curr_hit.query_name != next_hit.query_name \
           or curr_hit.target_name != next_hit.target_name \
           or not (0 <= query_hits_dist <= max_hits_dist) \
           or not (0 <= target_hits_dist <= max_hits_dist):
            united_hits.append(curr_hit)
            curr_hit = processed_hits[i]
        else:
            curr_hit.query_end = next_hit.query_end
            curr_hit.target_end = next_hit.target_end

    united_hits.append(curr_hit)
    return united_hits


def _complement_hit(hit):
    query_name = hit.query_name[:-1]
    query_name += ["+", "-"][hit.query_name.endswith("+")]
    query_start = hit.query_len - hit.query_end
    query_end = hit.query_len - hit.query_start

    target_name = hit.target_name[:-1]
    target_name += ["+", "-"][hit.target_name.endswith("+")]
    target_start = hit.target_len - hit.target_end
    target_end = hit.target_len - hit.target_start

    complement_hit = PafHit(
        query_name, hit.query_len, query_start, query_end,
        target_name, hit.target_len, target_start, target_end)

    complement_hit.id = -hit.id
    return complement_hit



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
            self.query_name, self.query_len,
            self.query_start, self.query_end)
        target_info = "{}\t{}\t{}\t{}".format(
            self.target_name, self.target_len,
            self.target_start, self.target_end)

        return query_info + target_info


def process_raw_hits(raw_hits, min_hit_length=20000):
    processed_hits = []

    for raw_hit in raw_hits:
        if raw_hit.query_hit_length() < min_hit_length:
            continue

        if raw_hit.target_hit_length() < min_hit_length:
            continue

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
    processed_hits.sort(key=lambda hit: (-hit.query_len,
                                         hit.query_name,
                                         hit.query_start))

    united_hits = []
    curr_hit = processed_hits[0]
    i = 1

    while i < len(processed_hits):
        next_hit = processed_hits[i]

        query_hits_dist = next_hit.query_start - curr_hit.query_end
        target_hits_dist = next_hit.target_start - curr_hit.target_end
        dist_diff = abs(query_hits_dist - target_hits_dist)

        if curr_hit.query_name != next_hit.query_name \
           or curr_hit.target_name != next_hit.target_name \
           or query_hits_dist < 0 or target_hits_dist < 0 \
           or dist_diff > 1000:
            united_hits.append(curr_hit)
            curr_hit = processed_hits[i]
        else:
            curr_hit.query_end = next_hit.query_end
            curr_hit.target_end = next_hit.target_end
        i += 1

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

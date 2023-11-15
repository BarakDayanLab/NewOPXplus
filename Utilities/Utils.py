import math
import json
import numpy as np
import itertools


class Utils:

    def __init__(self):
        raise Exception('This is a static class. No need to instantiate')
        pass

    @staticmethod
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    @staticmethod
    def moving_average(a, n=3):
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n

    @staticmethod
    def gaussian(x, mu, sigma):
        return 1 / sigma / np.sqrt(2 * np.pi) * np.exp(-1.0 / 2.0 * ((x - mu) / sigma) ** 2)

    @staticmethod
    def gaussian2(x, mu, sigma):
        return 0.2 * np.exp(-1.0 / 2.0 * ((x - mu) / sigma) ** 2)

    @staticmethod
    def gauss_adaptive(amplitude, length):
        t = np.linspace(-length / 2, length / 2, length)
        sigma = length / 6  # defining the width of the pulse to match the length desired
        gauss_wave = amplitude * np.exp(-(t ** 2) / (2 * sigma ** 2))
        return [float(x) for x in gauss_wave]

    @staticmethod
    def calc_cmat(correction_vars):
        """
        calculating the correction matrix required for IQ mixer using the variable \theta   k
        :param correction_vars: array of correction variables
        theta = correction_vars[0]
        theta = correction_vars[1]
        :return: correction matrix
        """
        theta = correction_vars[0]
        k = correction_vars[1]
        R = [[np.sin(theta), np.cos(theta)], [np.cos(theta), np.sin(theta)]]
        c = [[k, 0], [0, 1 / k]]
        return np.matmul(c, R).flatten().tolist()

    @staticmethod
    def most_common(self, lst):
        # get an iterable of (item, iterable) pairs
        SL = sorted((x, i) for i, x in enumerate(lst))
        # print 'SL:', SL
        groups = itertools.groupby(SL, key=operator.itemgetter(0))

        # auxiliary function to get "quality" for an item
        def _auxfun(g):
            item, iterable = g
            count = 0
            min_index = len(lst)
            for _, where in iterable:
                count += 1
                min_index = min(min_index, where)
            if count > 20:
                self.logger.debug('item %r, count %r, minind %r' % (item, count, min_index))
            return count, -min_index

        # pick the highest-count/earliest item
        return max(groups, key=_auxfun)[0]

    @staticmethod
    def amplitudeMultiplierToDBm(ratio):
        return 10 * np.log10(ratio)

    @staticmethod
    def dbmToP(dbm):
        return 10**(dbm/10)  # return mW

    @staticmethod
    def merge_jsons(json1, json2):
        """
        Merges json2 onto json1.
        - If keys are in json2 and missing in json1, they will be added
        - If keys exist in json1 and exist in json2, they will be overriden
        :param json1: Target Json
        :param json2: Source Json
        :return: the merged dictionary

        NOTE:

        The difference between JSON and Python dictionary:
        - The keys in JSON can be only strings.
        - The keys in the dictionary can be any hashable object.
        - In JSON, the keys are sequentially ordered and can be repeated.
        - In the dictionary, the keys cannot be repeated and must be distinct
        """

        obj1 = json.dumps(json1)
        obj2 = json.dumps(json2)

        dict1 = json.loads(obj1)
        dict2 = json.loads(obj2)

        merged_dict = {**dict1, **dict2}
        return merged_dict

    @staticmethod
    def merge_multiple_jsons(array_of_jsons):
        if len(array_of_jsons)==1:
            return array_of_jsons[0]

        # Iterate over all jsons in array and merge each one on top of the prev merger:
        result_json = array_of_jsons[0]
        for json_obj in array_of_jsons[1:]:
            result_json = Utils.merge_jsons(result_json, json_obj)
        return result_json

    @staticmethod
    def verify_keys_for_case_sensitivity(array_of_dictionaries):
        """
        Method for running sanity on keys in dictionaries - should not have the same keys with different case-sensitivity
        :param array_of_dictionaries: array of dictionaries
        :return: N/A

        i.e. it will generate a warning if we have 'MOT_Duration' in one dictionary and 'MOT_duration' in another.
        """

        # Concatenate all keys in all dictionaries
        keys = []
        for dic in array_of_dictionaries:
            keys = keys + list(dic.keys())

        # Iterate over all keys combinations and cross-check each pair
        for pair in itertools.combinations(keys, 2):
            if pair[0].lower() == pair[1].lower() and pair[0] != pair[1]:
                print('\033[91m' + f'Warning: {pair[0]} is not the same as {pair[1]}' + '\033[0m')
        pass

    @staticmethod
    def bucket_timetags_SIMPLE(timetags, window_size, buckets_number=1, filter=None, start_time=-math.inf, end_time=math.inf):
        """
        This function takes a sequence of time-tags and divides thjem into buckets.
        It does so scanning all time-tags from start to end and placing each timetag into the relevant bin.
        Bins are opened dynamically (e.g. if all time tags fall into the same bin, only one will return.
        The function returns both the counts and the time-tags themselves.

        :param timetags: An array of time tags we need to put into buckets
        :param window_size: The window size for the binning action
        :param buckets_number: Pre-allocated the buckets. If no time-tags for a specific bucket, it will return 0
                               If this param is None, there will be dynamic allocation until the last relevant timetag
        :param filter: A filter to apply for each time-tag - whether to count it or not, or weight - can be "bucket index"
        :param start_time: Filter out time-tags that come before this start time
        :param end_time: Filter out time-tags that come after this end time

        :return: Two arrays - (a) the buckets counts AND (b) the buckets with timetags
        """

        # TODO: a) Set local start/end times
        # TODO: b) better name for abs/loc timings - "folded" ?
        # TODO: c) add to buckets the folded or absolute?
        # TODO: d) different bucketing conditions? to different bucket results?

        # Initialize the data structures
        buckets_counts = [0] * buckets_number
        buckets_content = [[]]
        for x in range(0, buckets_number-1):
            buckets_content.append([])

        # Perform the bucketing
        for time_tag in timetags:
            if time_tag > start_time and time_tag < end_time:
                time_tag_in_sequence = time_tag % window_size
                seq_num = (time_tag - 1) // window_size
                increment = 1 if filter is None else filter[time_tag_in_sequence]

                # If we skipped few buckets (no time tags for them, create zero-count/zero-content buckets
                for i in range(len(buckets_counts), seq_num+1):
                    buckets_counts.append(0)
                    buckets_content.append([])
                buckets_counts[seq_num] += increment
                buckets_content[seq_num].append(time_tag)

        return buckets_counts, buckets_content

    @staticmethod
    def bucket_timetags(timetags, window_size, buckets_number=1, filter=None, start_time=-math.inf, end_time=math.inf, filters=None):
        """
        This function takes a sequence of time-tags and divides thjem into buckets.
        It does so scanning all time-tags from start to end and placing each timetag into the relevant bin.
        Bins are opened dynamically (e.g. if all time tags fall into the same bin, only one will return.
        The function returns both the counts and the time-tags themselves.

        :param timetags: An array of time tags we need to put into buckets
        :param window_size: The window size for the binning action
        :param buckets_number: Pre-allocated the buckets. If no time-tags for a specific bucket, it will return 0
                               If this param is None, there will be dynamic allocation until the last relevant timetag
        :param filters: A filter to apply on each time-tag - (a) whether to count it or not, (b) weight
        :param start_time: Filter out time-tags that come before this start time
        :param end_time: Filter out time-tags that come after this end time

        :return: Two arrays - (a) the buckets counts AND (b) the buckets with timetags
        """

        # Initialize the data structures
        buckets = {}
        for s in filters:
            if len(s["filter"]) < window_size:
                # TODO: Pad it or raise an exception? Add "strict" mode
                #raise Exception('Filter size must be equal or greater than window size!')
                gap = window_size - len(s["filter"])
                s["filter"] = np.pad(s["filter"], (0,gap))
            if len(s["filter"]) > window_size:
                print(f'Warning: filter size {len(s["filter"])} is larger than window size ({window_size})')

            #counts = [0] * buckets_number
            counts = np.zeros(buckets_number)
            content = [[]]
            for x in range(0, buckets_number-1):
                content.append([])
            name = s["name"]
            buckets[name] = {"name": name, "counts": counts, "content": content}

        # Perform the bucketing
        for absolute_time_tag in timetags:
            if absolute_time_tag > start_time and absolute_time_tag < end_time:
                relative_time_tag = absolute_time_tag % window_size
                seq_num = (absolute_time_tag - 1) // window_size

                # Decide with what bucket-set we are working
                for s in filters:
                    inc = s["filter"][relative_time_tag]
                    if inc == 0:
                        continue
                    b = buckets[s["name"]]

                    # If we skipped few buckets (no time tags for them, create zero-count/zero-content buckets
                    for i in range(len(b["counts"]), seq_num+1):
                        #b["counts"].append(0)
                        b["counts"] = np.append(b["counts"], 0)
                        b["content"].append([])
                    b["counts"][seq_num] += inc
                    b["content"][seq_num].append(absolute_time_tag)

        return buckets

    @staticmethod
    def bucket_test():

        time_tags = [1, 4, 15, 22, 24, 98, 99, 100]
        # s_filter = [i//50 for i in range(0, 100)]
        # n_filter = [1-(i//50) for i in range(0, 100)]
        s_filter = [0,1,0,1,0,1,0,1,0,1]  # Bucket of only odd timetags
        n_filter = [1,1,1,1,1,1,1,1,1,1]  # Bucket of all time tags

        buckets = Utils.bucket_timetags(
            timetags=time_tags,
            window_size=10,
            buckets_number=6,
            filter=None,
            start_time=1,
            end_time=100,
            filters=[{"name": "S", "filter": s_filter}, {"name": "N", "filter": n_filter}])

        pass

import os
import re
import math
import json
import uuid
import cv2, glob
import numpy as np
from scipy import ndimage
import itertools
from pkgutil import iter_modules
import pathlib
import collections.abc
import inspect
import subprocess
import matplotlib as mpl
import matplotlib.pyplot as plt

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
    def running_average(values, new_value, counter):
        return (values * (counter - 1) + new_value) / counter

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
    def solve_quadratic_equation(a, b, c):

        d = b ** 2 - 4 * a * c  # discriminant

        if d < 0:
            return []
        elif d == 0:
            x = (-b + math.sqrt(b ** 2 - 4 * a * c)) / 2 * a
            return np.array(x)
        else:
            x1 = (-b + math.sqrt((b ** 2) - (4 * (a * c)))) / (2 * a)
            x2 = (-b - math.sqrt((b ** 2) - (4 * (a * c)))) / (2 * a)
            return np.array([x1, x2])

    @staticmethod
    def calc_cmat(correction_vars):
        """
        Calculating the correction matrix required for IQ mixer using the variable \theta   k

        :param correction_vars: array of correction variables
        theta = correction_vars[0]
        theta = correction_vars[1]
        :return: correction matrix
        """
        theta = correction_vars[0]
        k = correction_vars[1]
        r = [[np.sin(theta), np.cos(theta)], [np.cos(theta), np.sin(theta)]]
        c = [[k, 0], [0, 1 / k]]
        return np.matmul(c, r).flatten().tolist()

    @staticmethod
    def most_common(self, lst):
        # get an iterable of (item, iterable) pairs
        sl = sorted((x, i) for i, x in enumerate(lst))
        # print 'SL:', SL
        groups = itertools.groupby(sl, key=operator.itemgetter(0))

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
    def data_to_list(data):
        """
        If data is ndarray --> list
        If data is int/float --> list
        If data is list --> list
        """
        if type(data) is np.ndarray:
            return data.tolist()
        elif type(data) != list and (type(data) == int or type(data) == float):
            return [data]
        return data

    @staticmethod
    def merge_jsons(json1, json2):
        """
        Merges json2 onto json1.
        - If keys are in json2 and missing in json1, they will be added
        - If keys exist in json1 and exist in json2, they will be overridden
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
        """
        Given an array of jsons, we merge them one onto the other, going from left to right
        This means that the first json in the array is the "weakest" and the last is the "strongest"
        """
        if len(array_of_jsons) == 1:
            return array_of_jsons[0]

        # Iterate over all jsons in array and merge each one on top of the prev merger:
        result_json = array_of_jsons[0]
        for json_obj in array_of_jsons[1:]:
            result_json = Utils.merge_jsons(result_json, json_obj)
        return result_json

    @staticmethod
    def _is_array(obj):
        """
        Checks if a given object is an array or np.array (not a scalar, float or string...)
        :param obj: Object to check
        :return: True/False
        """
        return not isinstance(obj, str) and isinstance(obj, (collections.abc.Sequence, np.ndarray))

    @staticmethod
    def verify_keys_for_case_sensitivity(array_of_dictionaries):
        """
        Method for running sanity on keys in dictionaries - shouldn't have the same keys with different case-sensitivity
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
    def bucket_timetags_SIMPLE(time_tags, window_size, buckets_number=1, filter=None, start_time=-math.inf, end_time=math.inf):
        """
        This function takes a sequence of time-tags and divides thjem into buckets.
        It does so scanning all time-tags from start to end and placing each timetag into the relevant bin.
        Bins are opened dynamically (e.g. if all time tags fall into the same bin, only one will return.
        The function returns both the counts and the time-tags themselves.

        :param time_tags: An array of time tags we need to put into buckets
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
        for time_tag in time_tags:
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
    def append_fifo(arr, value, size):
        """
        Appends a new value to an array.
        If the array size becomes larger than size, we throw out the first element (e.g. FIFO)
        """
        new_arr = arr[-(size-1):] + [value]
        return new_arr

    @staticmethod
    def split_into_value_ranges(values, value_ranges, default=False):
        """
        Split a list of values into the given ranges.
        Returns arrays with relevant values.
        All values that do not fall into any of the ranges go into a default bucket

        TODO: implement this more efficiently, no append, use fact it's sorted, sort ranges
        """
        # Initialize results array of arrays
        results = []
        for i in range(0, len(value_ranges)):
            results.append([])

        # Start going over values
        values.sort()
        for value in values:
            for i, value_range in enumerate(value_ranges):
                if value > value_range[0] and value < value_range[1]:
                    results[i].append(value)

        return results

        pass

    @staticmethod
    def elements_in_even_places(values):
        """
        Returns all elements which are in even indices. First element is in place 0, then 2, etc.
        If values is an empty array, and empty array is returned
        """
        return list(values[::2])

    @staticmethod
    def elements_in_odd_places(values):
        """
        Returns all elements which are in even indices. First element is in place 0, then 2, etc.
        If values is an empty array, and empty array is returned
        """
        return list(values[1::2])

    @staticmethod
    def split_to_even_odd_elements(values):
        """
        Return two arrays: first - the even-indexed elements & second - the odd-indexed elements
        """
        even_indexed_elements = Utils.elements_in_even_places(values)
        odd_indexed_elements = Utils.elements_in_odd_places(values)
        return even_indexed_elements, odd_indexed_elements

    @staticmethod
    def tail_contains_non_zeros(values, percent):
        """
        Checks that the tailing end of the array has some values that are non-zero
        TODO: add an option to have "head"/"trail"
        """
        percent = percent/100 if percent > 100 else percent
        index = int((1 - percent) * len(values))
        tail_array = np.array(values[index:])
        non_zero = sum(tail_array) != 0
        return non_zero

    @staticmethod
    def bin_values(values, bin_size, num_of_bins):
        """
        Create a histogram of time-tags with <Num-of-bins> bins and each bin holds values from 0-<bin_size-1>
        """
        binned_time_tags = np.zeros(num_of_bins)
        for x in values:
            binned_time_tags[(x - 1) // bin_size] += 1
        return binned_time_tags

    @staticmethod
    def bucket_timetags(timetags, window_size, buckets_number=1, filter=None, start_time=-math.inf, end_time=math.inf, filters=None):
        """
        This function takes a sequence of time-tags and divides them into buckets.
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
                # raise Exception('Filter size must be equal or greater than window size!')
                gap = window_size - len(s["filter"])
                s["filter"] = np.pad(s["filter"], (0, gap))
            if len(s["filter"]) > window_size:
                print(f'Warning: filter size {len(s["filter"])} is larger than window size ({window_size})')

            # counts = [0] * buckets_number
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
                        # b["counts"].append(0)
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
        s_filter = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]  # Bucket of only odd timetags
        n_filter = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # Bucket of all time tags

        buckets = Utils.bucket_timetags(
            timetags=time_tags,
            window_size=10,
            buckets_number=6,
            filter=None,
            start_time=1,
            end_time=100,
            filters=[{"name": "S", "filter": s_filter}, {"name": "N", "filter": n_filter}])

        pass

    # ---------------------------------------------------
    # Python Code Inspection related utilities
    # ---------------------------------------------------

    @staticmethod
    def type_string_to_type_function(type_str):
        if type_str == 'int' or type_str == 'integer':
            return int
        elif type_str == 'float':
            return float
        elif type_str == 'str' or type_str == 'string':
            return str
        else:
            return None  # This will probably cause the called to crash, but it's good - should check itself...
    @staticmethod
    def get_self_members(self_object, prefix=None, suffix=None, specific=None):
        """
        Return the members of the self object starting/ending with rule
        """
        _list = []
        for member in inspect.getmembers(self_object):
            # Ignore private/protected functions and methods that do not start with a underscore
            if not member[0].startswith('_') and not inspect.ismethod(member[1]):
                if prefix is not None and not member[0].startswith(prefix):
                    continue
                if suffix is not None and not member[0].endswith(suffix):
                    continue
                if specific is not None and member[0] not in specific:
                    continue
                _list.append(member[0])
        return _list

    @staticmethod
    def get_self_members_and_values(self_object, prefix=None, suffix=None, specific=None):
        """
        Return the members and values of the self object.
        It can be filtered using three options:
        - prefix - function will return only members that start with the prefix
        - suffix - function will return only members that end with the suffix
        - specific - an array of names - only these will return
        """
        _dict = {}
        for member in inspect.getmembers(self_object):
            # Ignore private/protected functions and methods that do not start with a underscore
            if not member[0].startswith('_') and not inspect.ismethod(member[1]):
                if prefix is not None and not member[0].startswith(prefix):
                    continue
                if suffix is not None and not member[0].endswith(suffix):
                    continue
                if specific is not None and member[0] not in specific:
                    continue
                _dict[member[0]] = member[1]
        return _dict

    @staticmethod
    def get_modules_in_package(package_dir=None):
        """
        Returns the names of the modules at a given package
        """
        if package_dir is None:
            package_dir = pathlib.Path(os.getcwd()).resolve()

        modules = []
        for (xxx, module_name, zzz) in iter_modules([package_dir]):
            modules.append(module_name)

        return modules

    # ---------------------------------------------------
    # Cavity related utilities
    # ---------------------------------------------------

    @staticmethod
    def cavity_transmission(f, k_ex, k_i, h):
        t = np.power(np.abs(1 + 2 * 1j * k_ex * (f - 1j * (k_i + k_ex)) /
                            (np.power((f - 1j * (k_i + k_ex)), 2) - np.power(h, 2))), 2)
        return t

    @staticmethod
    def cavity_reflection(f, k_ex, k_i, h):
        r = np.power(np.abs(2 * k_ex * h /
                            (np.power((1j * f + (k_i + k_ex)), 2) + np.power(h, 2))), 2)
        return r

    # ---------------------------------------------------
    # Network Related Utilities
    # ---------------------------------------------------

    @staticmethod
    def get_mac_address():
        node = uuid.getnode()
        arr = ['{:02x}'.format((node >> ele) & 0xff) for ele in range(0, 8 * 6, 8)][::-1]
        mac_str = ':'.join(arr)
        return mac_str

    # ---------------------------------------------------
    # Image/Video utilities
    # ---------------------------------------------------

    @staticmethod
    def create_mp4_from_path(input_path, output_path=None, file_name='video', fps=4):

        extension = 'mp4'  # 'mp4'

        if output_path is None:
            output_path = os.path.join(input_path, f'video.{extension}')

        img_array = []
        file_names = glob.glob(os.path.join(input_path, '*.bmp'))
        file_names.sort()
        the_size = (640, 480)
        for filename in file_names:
            img = cv2.imread(filename)
            height, width, layers = img.shape
            the_size = (width, height)
            img_array.append(img)

        fourcc = 'MP42'  # 'MP42', 'MP4V', 'DIVX'
        video_file_path = os.path.join(output_path, f'{file_name}.{extension}')
        out = cv2.VideoWriter(filename=video_file_path, fourcc=cv2.VideoWriter_fourcc(*fourcc), fps=fps, frameSize=the_size)
        for i in range(len(img_array)):
            out.write(img_array[i])

        out.release()
        pass

    @staticmethod
    def create_video_from_path(input_path, output_path=None, file_name='video'):
        """
        Creates a video from .bmp files in a given path
        - path - where the .bmp files are
        - save_file_path - Full path with name - where to save the generated video.
                           (if not provided, video is saved where the input images are as "video.avi")

        - Taken from: https://theailearner.com/2018/10/15/creating-video-from-images-using-opencv-python/
        """
        extension = 'avi'

        if output_path is None:
            output_path = input_path

        img_array = []
        file_names = glob.glob(os.path.join(input_path, '*.bmp'))
        file_names.sort()
        the_size = (640, 480)
        for filename in file_names:
            img = cv2.imread(filename)
            height, width, layers = img.shape
            the_size = (width, height)
            img_array.append(img)

        video_file_path = os.path.join(output_path, f'{file_name}.{extension}')
        fourcc = 'DIVX'
        out = cv2.VideoWriter(filename=video_file_path, fourcc=cv2.VideoWriter_fourcc(*fourcc), fps=3, frameSize=the_size)
        for i in range(len(img_array)):
            out.write(img_array[i])
        out.release()

    # ---------------------------------------------------
    # OS related utilities
    # ---------------------------------------------------

    def is_port_established(port):
        netstat_lines = Utils.get_netstat_lines()
        results = [line for line in netstat_lines if f':{port} ' in line and 'ESTABLISHED' in line]
        return len(results) > 0

    @staticmethod
    def get_netstat_lines():
        lines = []
        try:
            res = subprocess.run(['netstat', '-n'], stdout=subprocess.PIPE)
            lines = str(res.stdout).split('\\r\\n')
        except Exception as err:
            print(f'Failed to launch explorer at path. {err}')
        return lines

    @staticmethod
    def open_windows_explorer(folder):
        """
        Opens windows explorer at a given path
        """
        cmd = f'explorer /select,"{folder}"'
        try:
            subprocess.Popen(cmd)
        except Exception as err:
            print(f'Failed to launch explorer at path {folder}. {err}')
        pass

    @staticmethod
    def ensure_folder_exists(folder_path):
        """
        Creates a folder if it doesn't exist. Returns a flag indicating if it existed before
        :param folder_path: folder path to create/check
        :return: The existence of the file before the action
        """
        existed = True
        if not os.path.exists(folder_path):
            os.makedirs(folder_path, exist_ok=True)
            existed = False
        return existed

    @staticmethod
    def remove_empty_experiment_folders(root_folder):
        """
        Cleans-up empty experiment folders. These are folders that contain only a single file: log.txt
        """

        folders = Utils.get_folders_in_path(root_folder)

        # Iterate over all folders and check which one has only a single file: "log.txt"
        folders_to_remove = []
        for folder in folders:
            folder = os.path.join(root_folder, folder)
            all_files = Utils.get_files_in_path(folder, return_full_path=False)
            if len(all_files) == 1 and all_files[0] == 'log.txt':
                folders_to_remove.append(folder)

        # Remove the folders with only "log.txt" file
        for folder in folders_to_remove:
            try:
                os.remove(os.path.join(folder, 'log.txt'))
                os.rmdir(folder)
                print(f'- Removed {folder}')
            except Exception as err:
                print(f'*** Failed to remove folder {folder}. {err}')

        pass

    @staticmethod
    def get_files_in_path(path, opt_in_filter=None, opt_out_filter=None, return_full_path=False):
        """
        Returns all files that are in a given path.
        - opt_out_filter - Exclude files that have the given token
        - opt_in_filter - Include only files that have the given token
        - return_full_path - full-path/file-name or only file-name
        """
        files = []
        for f in os.listdir(path):
            if os.path.isfile(os.path.join(path, f)):
                if opt_out_filter is not None and opt_out_filter in f:
                    continue
                if opt_in_filter is None or (opt_in_filter is not None and opt_in_filter in f):
                    if return_full_path:
                        files.append(os.path.join(path, f))
                    else:
                        files.append(f)
        return files

    @staticmethod
    def get_folders_in_path(path):
        folders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
        return folders

    @staticmethod
    def load_json_from_file(file_full_path):
        try:
            f = open(file_full_path)
            json_obj = json.load(f)
            f.close()
            return json_obj
        except Exception as err:
            pass
        return None

    # ---------------------------------------------------
    # Other Utilities
    # ---------------------------------------------------

    @staticmethod
    def generate_UUID():
        myuuid = uuid.uuid4()
        return myuuid

    @staticmethod
    def box_string(msg):
        rpt = '=*' * 40
        str = ('\n' + rpt + '\n')
        str += '*=\n'
        str += ('=*   ' + msg + '\n')
        str += '*=\n'
        str += (rpt + '\n')
        return str

    @staticmethod
    def _extract_key_value_pairs_from_file_name(file_name):
        # Remove extensions and change comma to semi-colon
        file_name = file_name.replace('.jpeg', '').replace('.bmp', '').replace(',', ';')

        key_and_values = file_name.split(';')
        key_value_pairs = []
        for key_str in key_and_values:
            value = key_str[key_str.find('=') + 1:]
            key = key_str[:key_str.find('=')]
            key_value_pair = (key, float(value))
            key_value_pairs.append(key_value_pair)
        return key_value_pairs

    @staticmethod
    def _extract_ppd_from_file_name(file_name):
        match = re.search(r'[-+]?\d*\.\d+|\d+', file_name)
        ppd = float(match.group())
        return ppd

    @staticmethod
    def _crop_around_center(self, image, center, crop_img_size):

        X_PIXEL_LEN = 1544
        Y_PIXEL_LEN = 2064

        X_UPPER_BOUND = int(center[0] + crop_img_size)
        X_LOWER_BOUND = int(center[0] - crop_img_size)
        if X_LOWER_BOUND < 0:
            X_LOWER_BOUND = 0
        if X_UPPER_BOUND > X_PIXEL_LEN:
            X_UPPER_BOUND = X_PIXEL_LEN

        Y_UPPER_BOUND = int(center[1] + crop_img_size)
        Y_LOWER_BOUND = int(center[1] - crop_img_size)
        if Y_LOWER_BOUND < 0:
            Y_LOWER_BOUND = 0
        if Y_UPPER_BOUND > Y_PIXEL_LEN:
            Y_UPPER_BOUND = Y_PIXEL_LEN

        # Crop the image
        img = image[Y_LOWER_BOUND:Y_UPPER_BOUND, X_LOWER_BOUND:X_UPPER_BOUND]
        return img

    @staticmethod
    def plot_cloud_position_over_time(images_folder, out_path, x_start, x_end, y_start, y_end,
                                      mm_to_pixel, start_image, skip_images, show, debug):

        mpl.rcParams['figure.dpi'] = 300

        MM_TO_PXL = 10/1000  # 10 / 970

        # Load the background image
        backgroundImg = cv2.imread(os.path.join(images_folder, 'background.bmp'), 0)

        # Get the files we want to show
        files = [f for f in os.listdir(images_folder) if
                 os.path.isfile(os.path.join(images_folder, f)) and 'PrePulse_duration' in f]
        files = files[start_image::skip_images]

        fig, ax = plt.subplots(nrows=1, ncols=1)

        all_images = []
        ppds = []
        y_centers = []
        for index, f in enumerate(files):

            # Extract the ppd value
            ppd = Utils._extract_ppd_from_file_name(f)
            ppds.append(ppd)

            # Load the image and subtract the background from it
            full_file_name = os.path.join(images_folder, f)
            ImgToFit = cv2.imread(full_file_name, 0)
            ImgToFit = cv2.subtract(ImgToFit, backgroundImg)

            # Crop the image
            ImgToFit = ImgToFit[y_start:y_end, x_start:x_end]

            # Append the cloud center (found by max intensity)
            img_y_max_index = np.argmax(np.sum(ImgToFit, axis=1))
            y_centers.append(img_y_max_index)

            # Find center of mass
            # res = ndimage.measurements.center_of_mass(ImgToFit)
            # ImgToFit = cv2.circle(ImgToFit, (int(res[0]), int(res[1])), radius=3, color=(0, 0, 255), thickness=-1)

            # Add PrePulseDuration value on image
            image_with_text = cv2.putText(ImgToFit, text=f'TOF={ppd} ms', org=(50,50), fontFace=cv2.FONT_HERSHEY_SIMPLEX,
                                fontScale=1, color=(255, 0, 0), thickness=2, lineType=cv2.LINE_AA)

            # Append image to array of images
            all_images.append(image_with_text)

            if debug:
                cv2.imshow(f'ppd={ppd}, y={img_y_max_index}', ImgToFit)
                cv2.waitKey(0)
            pass

        # Concatenate all images horizontally
        concatenated_images = cv2.hconcat(all_images)

        delta = y_centers[0]
        y_top = 0 + delta
        y_bottom = -(y_end - y_start) + delta

        # Translate pixels to millimeters
        y_top *= mm_to_pixel
        y_bottom *= mm_to_pixel

        ax.imshow(concatenated_images, extent=[0, 40, y_bottom, y_top], aspect=1)  # Min/Max, Max/Min

        # Graph Cosmetics
        plt.title('Cloud positions in different timings', fontsize=6)
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.ylabel('Z Position (mm)', fontsize=6)
        plt.tick_params(axis='y', labelsize=6)

        # Set y-axis ticks with 0.5 intervals and Celsius labels
        yticks = np.arange(4, -2, -1)
        plt.yticks(yticks)
        #plt.yticks(yticks, labels=[f'{(5 / 9) * (f - 32):.1f}°C' for f in yticks])

        if out_path is not None:
            plt.savefig(out_path, dpi=300, bbox_inches=None, pad_inches=0.05)
            # plt.savefig(out_path, dpi=300, bbox_inches='tight', pad_inches=0.2)
            # plt.savefig(out_path, dpi=300)

        if show:
            plt.show()

        pass


if __name__ == "__main__":

    #images_folder = r'C:\temp\bitmaps_for_sunrise_chart'
    #images_folder = r'C:\temp\bitmaps_for_sunrise_chart_2'
    images_folder = r'C:\temp\bitmaps_for_sunrise_20240903_181633'

    out_image = os.path.join(images_folder, 'MOT_trajectory.jpg')
    Utils.plot_cloud_position_over_time(images_folder=images_folder, out_path=out_image,
                # x_start=800, x_end=1200, y_start=1000, y_end=1500,
                x_start=800, x_end=1200, y_start=600, y_end=1200,
                mm_to_pixel=10/970, start_image=0, skip_images=2,
                show=True, debug=False)

    # Test
    # Utils.remove_empty_experiment_folders(r'C:\temp\refactor_debug\Experiment_results\Temperature')
    pass
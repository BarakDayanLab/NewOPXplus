import os
import time
import struct
import numpy as np
from io import BytesIO
from Utilities.Utils import Utils

# -------------------------------------------------------------------------------------------------------
# TODO:
# -------------------------------------------------------------------------------------------------------
# 1. Write the file in a more compact way - the len numbers can be 4 bytes and not 8 (int and not float)
# 2. In the detectors data, we can write integers instead of doubles. Not so in the Flr data
# 3. The read process can also probably be done faster - instead of reading each value and appending it
# 4. Write timestamps at the beginning of each save
# 5. Add "Seek" method - to seek to a specific time stamp
# 6. Make this class a generic one for all stream purposes - fetching from OPX, etc.
# 7. At the end of experiment - copy it from local to experiment folder
# -------------------------------------------------------------------------------------------------------


class BDStreams:

    def __init__(self, streams=None, save_path=None, save_streams=False):
        self.save_path = save_path
        self.save_name = os.path.join(self.save_path, 'raw_streams.dat')
        self.save_streams = save_streams

        self.streams_defs = streams
        self.number_of_rows_saved = 0

        pass

    def __del__(self):
        pass

    def set_streams_definitions(self, streams):
        self.streams_defs = streams
        pass


    def normalize_stream(self, stream_data, detector_index, detector_delays, M_window):
        """
        Given a stream data, normalize it:
        - Remove the first element - the count
        - Remove time-tags who are exactly on the M_window
        - Remove junk values coming from the OPX (9999984)
        - Remove time-tags that after adding the delay will be outside the window
        """
        stream_len = stream_data[0]
        delay = detector_delays[detector_index]

        # Iterate starting from the first position after the count and build the normalized stream
        normalized = []
        for tt_index in range(1, stream_len+1):
            data = stream_data[tt_index]
            if data % M_window == 0 or data == 9999984:
                continue

            delayed_time_tag = data + delay
            if delayed_time_tag < M_window:
                normalized.append(delayed_time_tag)

        normalized.sort()
        return normalized

    def clean_streams(self):
        for stream in self.streams_defs.values():
            stream['all_rows'] = []
            stream['results'] = []

    # TODO: complete this:
    def load_entire_folder(self, playback_files_path):
        """
        Iterates over all files in folder and loads them
        """

        # Clean all streams
        self.clean_streams()

        # Get all files in playback folder
        playback_files = [f for f in os.listdir(playback_files_path) if os.path.isfile(os.path.join(playback_files_path, f))]

        # Iterate over all files in folder
        for playback_file in playback_files:
            self.load_streams_enhanced(os.path.join(playback_files_path, playback_file))


        pass

    def load_streams_enhanced(self, data_file):

        # Open the binary file
        with open(data_file, 'rb') as file:

            try:
                bytes = file.read(4)
                bytes_unpacked = struct.unpack('>i', bytes)[0]
                current_time = int(bytes_unpacked)

                number_of_streams = int(struct.unpack('>b', file.read(1))[0])
                for i in range(0, number_of_streams):

                    # Get stream name len and then name
                    name_len = struct.unpack(f'>I', file.read(4))[0]
                    name = struct.unpack(f'>{name_len}s', file.read(name_len))[0]
                    name = name.decode()

                    # Get data len and then data
                    data_len = struct.unpack('>I', file.read(4))[0]
                    data_bytes = file.read(data_len)

                    # Reconstruct the numpy array from bytes
                    load_bytes = BytesIO(data_bytes)
                    loaded_np = np.load(load_bytes, allow_pickle=True)

                    # Append the data we got to all rows
                    stream = self.streams_defs[name]
                    if 'all_rows' not in stream:  # If it's the first results we're adding, create a new array
                        stream['all_rows'] = [loaded_np]
                    else:
                        stream['all_rows'].append(loaded_np)

                    # TODO: we should append these - to create the sequence of time stamps recorded. Now we override.
                    stream['timestamp'] = current_time

            except Exception as err:
                print(err)
                pass
        pass

    def load_streams(self, data_file):

        # Create an array of stream names (because we'll be working with indices)
        #stream_names = list(self.streams_defs.keys())

        # Get only those streams that their configuration indicates they need to be saved
        streams_to_load = [s for s in self.streams_defs.values() if 'save_raw' in s and s['save_raw']]

        # Open the binary file
        with open(data_file, 'rb') as file:

            try:
                bytes = file.read(4)
                bytes_unpacked = struct.unpack('>i', bytes)[0]
                current_time = int(bytes_unpacked)
                number_of_streams = int(struct.unpack('>b', file.read(1))[0])
                for i in range(0, number_of_streams):
                    if i == 8:
                        stream_name = 'FLR_measure'
                    else:
                        stream_name = f'Detector_{i+1}_Timetags'  # TODO: need to do something smarter in the save - and then here

                    stream = self.streams_defs[stream_name]
                    #stream = self.streams_defs[stream_names[i]]
                    stream_data_len = struct.unpack('>H', file.read(2))[0]
                    type_func = Utils.type_string_to_type_function(stream['type'])

                    # Read the data of stream
                    results = []
                    for j in range(0, stream_data_len):
                        # Get the value
                        sz = struct.calcsize(f'>{stream["binary"]}')
                        val = struct.unpack(f'>{stream["binary"]}', file.read(sz))[0]
                        # Cast it
                        val = type_func(val)
                        # Keep it
                        results.append(val)

                        # If it's the first results we're adding, create a new array
                        if 'all_rows' not in stream:
                            stream['all_rows'] = []

                    # Append it to all rows
                    stream['all_rows'].append(results)
                    stream['timestamp'] = current_time
            except Exception as err:
                print(err)
                pass
        pass

    def _data_to_list(self, data):
        """
        If data is ndarray --> list
        If data is int --> list
        If data is list --> list
        """
        if type(data) is np.ndarray:
            return data.tolist()
        elif type(data) != list and (type(data) == int or type(data) == float):
            return [data]
        return data

    def save_streams_enhanced(self):
        """

                        +-----------------+-------------+
        Stream Name ==> | Stream Name Len | Stream Name |
                        +-----------------+-------------+

                        +-----------------+-------------+
        Stream Data ==> | Stream Data Len | Stream Data +
                        +-----------------+-------------+

        +-----------+-------------+---------------+---------------+-------+---------------+---------------+
        + Timestamp | Num Streams | Stream-1 Name | Stream-1 Data | ..... | Stream-N Name | Stream-N Data |
        +-----------+-------------+---------------+---------------+-------+---------------+---------------+

        """
        if not self.save_streams:
            return

        if self.streams_defs is None:
            print('No streams defined - nothing to save! Ignoring')
            return

        current_time = time.time()
        bytes_array = bytearray()

        try:
            # Get only those streams that their configuration indicates they need to be saved
            streams_to_save = [s for s in self.streams_defs.values() if 'save_raw' in s and s['save_raw']]

            # Start packing the data with 'b' (1-byte=unsigned-char) for number of streams
            bytes_array += struct.pack('>ib', int(current_time), len(streams_to_save))

            # Iterate over all streams and pack their name and data
            for stream in streams_to_save:

                # Pack the name (I=unsigned int for len and then name)
                name = stream['name']
                bytes_array += struct.pack(f'>I{len(name)}s', len(name), bytes(name, 'utf-8'))

                # Pack the data
                data = [] if stream['handler'] is None else stream['results']

                # Save in to BytesIo buffer
                np_bytes = BytesIO()
                np.save(np_bytes, data, allow_pickle=True)

                # get bytes value
                np_bytes = np_bytes.getvalue()
                bytes_array += (struct.pack(f'>I', len(np_bytes)) + np_bytes)

        except Exception as err:
            print(f'Failed to save raw data: {err}')

        time_formatted = time.strftime("%Y%m%d_%H%M%S")
        save_name = os.path.join(self.save_path, f'{time_formatted}_streams.dat')

        with open(save_name, "wb") as file:
            file.write(bytes_array)

        self.number_of_rows_saved += 1

        total_prep_time = time.time() - current_time
        pass

    def save_streams(self):
        """
        Iterate over all streams and save to a file with time-stamp
        Structure of binary file:
        - Records, each holding the Number of Streams
        - Then, for each stream, the size of it
        - The stream data

        The reasonable time for a single file save is 1-2 ms.
        +-----------+-------------+--------------+---------------+-------+--------------+---------------+
        + Timestamp + Num Streams | Stream-1 Len | Stream 1 Data | ..... | Stream-N Len | Stream-N Data |
        +-----------+-------------+--------------+---------------+-------+--------------+---------------+

        """
        if self.streams_defs is None:
            print('No streams defined - nothing to save! Ignoring')
            return

        current_time = time.time()

        try:
            # Get only those streams that their configuration indicates they need to be saved
            streams_to_save = [s for s in self.streams_defs.values() if 'save_raw' in s and s['save_raw']]

            # Start packing the data with 'B'=1-byte for number of streams
            bytes_array = struct.pack('>ib', int(current_time), len(streams_to_save))
            for stream in streams_to_save:

                data = [] if stream['handler'] is None else stream['results']

                # TODO: when we save the data correctly, we should not need this, as everything will be arrays
                # TODO: unlike now where FLR is saved as a single float value
                data = self._data_to_list(data)

                # Pack the data with: (a) 'H'=2-bytes-Unsigned Short for stream size (b) 'I'=4-bytes-Unsigned int for stream data
                data = [len(data)] + data
                data_packed = struct.pack(f'>H{len(data) - 1}{stream["binary"]}', *data)
                bytes_array = bytes_array + data_packed

        except Exception as err:
            print(f'Failed to save raw data: {err}')

        time_formatted = time.strftime("%Y%m%d_%H%M%S")
        save_name = os.path.join(self.save_path, f'{time_formatted}_streams.dat')

        with open(save_name, "wb") as file:
            file.write(bytes_array)

        self.number_of_rows_saved += 1

        total_prep_time = time.time() - current_time
        pass

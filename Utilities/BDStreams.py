import os
import time
import struct
import numpy as np
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

    def __init__(self, streams=None, save_path=None):
        self.save_path = save_path
        self.save_name = os.path.join(self.save_path, 'raw_streams.dat')

        self.streams_defs = streams
        self.number_of_rows_saved = 0

        # The file handle for cases we would like to write all to a singel file
        self.file = None

        pass

    def __del__(self):
        self.close_file()

    def close_file(self):
        if self.file is not None:
            self.file.close()
            self.file = None

    def set_streams_definitions(self, streams):
        self.streams_defs = streams

    def clean_streams(self):
        for stream in self.streams_defs.values():
            stream['all_rows'] = []
            stream['results'] = []

    # TODO: complete this:
    def load_entire_folder(self, raw_streams_folder):

        # Clean all streams
        self.clean_streams()

        # Get all files in folder
        # TODO: ...

        # Iterate over files, and for each, run the following:
        self.load_streams(r'C:\temp\streams_raw_data_2\20231225_143721_streams.dat')

        pass

    def load_streams(self, data_file):

        # Create an array of stream names (because we'll be working with indices)
        stream_names = list(self.streams_defs.keys())

        # Open the binary file
        with open(data_file, 'rb') as file:

            number_of_streams = int(struct.unpack('>b', file.read(1))[0])
            for i in range(0, number_of_streams):
                stream = self.streams_defs[stream_names[i]]
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
        pass

    def load_streams_from_a_single_file(self):

        # Create an array of stream names (because we'll be working with indices)
        stream_names = list(self.streams_defs.keys())

        # Initialize the file position and size for the reading loop
        file_position = 0
        file_size = os.stat(self.save_name).st_size

        # Open the binary file
        with open(self.save_name, 'rb') as file:

            # Iterate while we still have something to read
            while file_position < file_size:

                number_of_streams = int(struct.unpack('>b', file.read(1))[0])
                for i in range(0, number_of_streams):
                    stream_name = stream_names[i]
                    stream = self.streams_defs[stream_name]
                    stream_data_len = int(struct.unpack('I', file.read(4))[0])
                    type_func = Utils.type_string_to_type_function(self.streams_defs[stream_name]['type'])

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

                    # If its the first results we're adding, create a new array
                    if 'all_rows' not in self.streams_defs[stream_name]:
                        self.streams_defs[stream_name]['all_rows'] = []

                    # Append it to all rows
                    self.streams_defs[stream_name]['all_rows'].append(self.streams_defs[stream_name]['all_rows'])

                file_position = file.tell()
        pass

    def save_streams(self):
        """
        Iterate over all streams and save to a file with time-stamp
        Structure of binary file:
        - Records, each holding the Number of Streams
        - Then, for each stream, the size of it
        - The stream data

        The reasonable time for a single file save is 1-2 ms.
        +-------------+--------------+---------------+-------+--------------+---------------+
        + Num Streams | Stream-1 Len | Stream 1 Data | ..... | Stream-N Len | Stream-N Data |
        +-------------+--------------+---------------+-------+--------------+---------------+
        """
        if self.streams_defs is None:
            print('No streams defined - nothing to save! Ignoring')
            return

        start_save_time = time.time()

        # Start packing the data with 'B'=1-byte for number of streams
        bytes_array = struct.pack('>b', len(self.streams_defs))
        all_data = []
        for stream in self.streams_defs.values():
            data = [] if stream['handler'] is None else stream['results']
            # TODO: when we save the data correctly, we should not need this, as everything will be arrays
            # TODO: unlike now where FLR is saved as a single float value
            if type(data) != list:
                data = [data]

            # Pack the data with: (a) 'H'=2-bytes-Unsigned Short for stream size (b) 'I'=4-bytes-Unsigned int for stream data
            data = [len(data)] + data
            data_packed = struct.pack(f'>H{len(data) - 1}{stream["binary"]}', *data)
            bytes_array = bytes_array + data_packed

            all_data = all_data + [len(data)] + data
        # Write number of streams in "line" of data and then the stream data itself
        all_data = [len(self.streams_defs)] + all_data

        time_formatted = time.strftime("%Y%m%d_%H%M%S")
        save_name = os.path.join(self.save_path, f'{time_formatted}_streams.dat')

        with open(save_name, "wb") as file:
            file.write(bytes_array)

        self.number_of_rows_saved += 1

        total_prep_time = time.time() - start_save_time

        pass

    # TODO: this method needs to be refactored - to use the code written in the save_stream() method
    def save_streams_append(self):
        """
        Iterate over all streams and append to file
        Structure of binary file:
        - Records, each holding the Number of Streams
        - Then, for each stream, the size of it
        - The sttream data
        """
        if self.streams_defs is None:
            print('No streams defined - nothing to save! Ignoring')
            return

        if self.file is None:
            self.file = open(self.save_name, 'ab')

        start_save_time = time.time()

        all_data = []
        for stream in self.streams_defs.values():
            data = [] if stream['handler'] is None else stream['results']
            # TODO: when we save the data correctly, we should not need this, as everything will be arrays
            # TODO: unlike now where FLR is saved as a single float value
            if type(data) != list:
                data = [data]
            all_data = all_data + [len(data)] + data
        # Write number of streams in "line" of data and then the stream data itself
        all_data = [len(self.streams_defs)] + all_data

        # TODO: make this more compact - like we do today in the save_file method - take struct.pack code from there.
        fmt = f'{len(all_data)}d'  # Double
        self.file.write(struct.pack(fmt, *all_data))

        # TODO: We can decide to flush every few cycles and not every save
        self.file.flush()

        self.number_of_rows_saved += 1

        total_prep_time = time.time() - start_save_time

        pass
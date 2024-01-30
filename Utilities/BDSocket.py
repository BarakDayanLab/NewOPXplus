import os
import json
import socket
import threading


class BDSocket:

    SERVER_IP = "127.0.0.1"  # server hostname or IP address
    PORT = 5050  # server port number

    def __init__(self, writeable):

        self.writeable = writeable

        # Load the connection definitions into map
        try:
            the_path = '.'
            the_file = os.path.join(the_path, 'connections.json')
            f = open(the_file)
            self.connections_map = json.load(f)
            f.close()
        except Exception as err:
            self._handle_error(f'Unable to open connections.json. Reason: {err}', True)
        pass

    def handle_client(self, client_socket, addr, connection_name, writeable):
        try:
            while True:
                # Receive client messages
                request = client_socket.recv(1024).decode("utf-8")

                # Get the relevant entry
                connection = self.connections_map[connection_name]

                if 'callback' in connection:
                    connection['callback'](request)
                    continue
                else:
                    # Assign the value received
                    writeable[connection['value_name']] = request

                if request.lower() == "close":
                    client_socket.send("closed".encode("utf-8"))
                    break
                print(f"Received: {request}")
                # convert and send accept response to the client
                response = "accepted"
                client_socket.send(response.encode("utf-8"))
        except Exception as e:
            print(f"Error when handling client: {e}")
        finally:
            client_socket.close()
            print(f"Connection to client ({addr[0]}:{addr[1]}) closed")

    def run_server(self):
        server_thread = threading.Thread(target=self._run_server)
        server_thread.start()

    # def _get_connection_by_key(self, key, value):
    #     """ Look for the connection that has the specific value in the key defined
    #         Return None if not found.
    #     """
    #     # Iterate over all connections
    #     for attribute, entry in self.connections_map.items():
    #         if value == entry[key] or value in entry[key]:
    #             return attribute, entry
    #     return None, None

    def _get_connection_entry(self, address, port):
        """ Look for the connection that has this address defined
            Return None if not found (e.g. - not an approved connection)
        """
        # Iterate over all connections
        for attribute, entry in self.connections_map.items():
            if '*' in entry['trusted_IP']:
                return entry
            if address in entry['trusted_IP']:
                return attribute, entry
        return None, None

    def _run_server(self):
        # Create a socket object
        server = None
        try:
            server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            # bind the socket to the host and port
            server.bind((BDSocket.SERVER_IP, BDSocket.PORT))
            # listen for incoming connections
            server.listen()
            print(f"Listening on {BDSocket.SERVER_IP}:{BDSocket.PORT}")

            while True:
                # accept a client connection
                client_socket, addr = server.accept()

                # Check if this is a connection we are ready to accept
                connection_name, connection_entry = self._get_connection_entry(addr[0], addr[1])
                if connection_entry is None:
                    print(f'Attempt to connect from {addr[0]}:{addr[1]}. Not approved!')
                    continue

                # Accept the connection & start the new thread to handle it
                print(f"Accepted connection from {addr[0]}:{addr[1]}")
                connection_entry['thread'] = threading.Thread(target=self.handle_client, args=(client_socket, addr, connection_name, self.writeable))
                connection_entry['thread'].start()

        except Exception as e:
            print(f"Error: {e}")
        finally:
            server.close()

    @staticmethod
    def run_client():
        """
        Run a test client - connect to the server and send a message
        """
        # Create a socket object
        client = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        # Establish connection with server
        client.connect((BDSocket.SERVER_IP, BDSocket.PORT))

        try:
            while True:
                # get input message from user and send it to the server
                msg = input("Enter message: ")
                client.send(msg.encode("utf-8")[:1024])

                # receive message from the server
                response = client.recv(1024)
                response = response.decode("utf-8")

                # If server sent us "closed" in the payload, we break out of
                # the loop and close our socket
                if response.lower() == "closed":
                    break

                print(f"Received: {response}")
        except Exception as e:
            print(f"Error: {e}")
        finally:
            # close client socket (connection to the server)
            client.close()
            print("Connection to server closed")

    def _handle_error(self, message, raise_exception=False):
        OK_COLOR = '\033[94m'
        ERR_COLOR = '\033[91m'
        cc = ERR_COLOR  # OK_COLOR
        message = f"{cc} NOTE: Error in connections: {message} {cc}"
        print(message)
        if raise_exception:
            raise Exception(f'Failed. Aborting.')


if __name__ == "__main__":

    # Run a test client
    BDSocket.run_client()

    pass

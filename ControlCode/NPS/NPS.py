"""This is a class to control the ePDU from Eaton.

Connection information as well as timeouts are saved at the file
    pointed to in the path variable below

If this file is changed, please make sure to put it in the python path
    for any environments that may use it ('/home/vfndev/anaconda3/envs/{env name}/lib/python3.7')    
"""

# python standard library
from telnetlib import Telnet
from configparser import ConfigParser

nps_config_file = "/home/vfndev/Documents/.NPS_config.ini"

class NPS:
    """A class and context manager for telnet control of a Pulizzi"""

    def __init__(self):
        """Constructor, takes the IP and the port of the Pulizzi

        Args:
            address   = the ip address of the Pulizzi
            port      = the telnet port of th Pulizzi
        """

        cp = ConfigParser()
        cp.read(nps_config_file)

        # define a timeout interval in seconds
        self.tmt = cp.getint("Settings", "timeout")

    def q_port(self, port):
        """Queries the on status of ports

        Args:
            port = can be an int representing a port, a list of ports, or 'all' for all ports
        Returns:
            dict = a dictionary with keys as port numbers, values as True/False for on/off
        """

        # if there's only one element in the list, convert it to an int
        if type(port) is list and len(port) == 1:
            try: port = int(port[0])
            except ValueError: raise ValueError("List elements should be a number 1 through 8.")

        # if port is an int, there's only port to check, so check it
        if type(port) is int:
            if port not in range(1, 9):
                raise ValueError("NPS port must be between 1 and 8")

            # send On/Off query request for given outlet (port)
            self.telnet.write("get PDU.OutletSystem.Outlet[{}].PresentStatus.SwitchOnOff\r\n".format(port).encode())

            # the two valid responses are 1 for on, 0 for off
            res = self.telnet.expect([bytes("1\r\r\npdu#0>", "utf-8"), bytes("0\r\r\npdu#0>", "utf-8")], self.tmt)
            if res[0] == 0: return {port:True}
            elif res[0] == 1: return {port:False}
            else: raise ConnectionError("Response to On/Off query not understood.")

        # in the 'all ports' case, make a list of the valid ports
        elif type(port) is str and port.lower() == "all":
            port = [1, 2, 3, 4, 5, 6, 7, 8]

        # if what remains isn't a list, throw an error
        elif type(port) is not list:
            raise ValueError("Please pass an int, a list of ints, or the string 'all'.")

        # a dictionary to hold the results
        res = {}

        # go through the ports, querying one by one
        for p in port:
            try: p = int(p)
            except ValueError: raise ValueError("Ports must be integers.")
            res[p] = self.q_port(p)[p]

        return res

    def s_port(self, requests):
        """Sets the given ports to on or off

        Args:
            requests = a dictionary where keys correspond to ports (or all), and values are True/False for On/Off
        Returns:
            None
        """

        # if all off or all on was sent, reformat requests
        if "all" in requests.keys():
            requests = {port : requests["all"] for port in [1, 2, 3, 4, 5, 6, 7, 8]}

        for port in requests:
            # check that port number is valid
            try:
                if int(port) not in range(1, 9):
                    raise ValueError("Ports must be between 1 and 8 (or 'all')")
            except ValueError: raise ValueError("Ports must be integers.")

            # check if request is on or off
            if requests[port]:
                self.telnet.write("set PDU.OutletSystem.Outlet[{}].DelayBeforeStartup 0\r\n".format(port).encode())
            else:
                self.telnet.write("set PDU.OutletSystem.Outlet[{}].DelayBeforeShutdown 0\r\n".format(port).encode())

            # we should get the pdu line start again
            res = self.telnet.expect([bytes("pdu#0>", "utf-8")], self.tmt)
            if res[0] == -1:
                msg = "Issue turning on port {}.".format(port)
                raise ConnectionError(msg)

    def __enter__(self):
        """Opens connection"""

        cp = ConfigParser()
        cp.read(nps_config_file)

        # Start telnet connection
        try: self.telnet = Telnet(cp.get("Connection Info", "address"),
            cp.getint("Connection Info", "port"), self.tmt)
        except: 
            raise ConnectionError("Connection failed.")

        # send login
        res = self.telnet.expect([bytes("Enter Login: ", 'utf-8')], self.tmt)
        if(res[0] == -1):
            raise ConnectionError("Connection timeout on startup. (No login prompt)")

        self.telnet.write("{}\r\n".format(cp.get("Connection Info", "login")).encode())

        # send password
        res = self.telnet.expect([bytes("Enter Password: ", 'utf-8')], self.tmt)
        if(res[0] == -1):
            raise ConnectionError("Connection timeout on startup. (No password prompt)")

        self.telnet.write("{}\r\n".format(cp.get("Connection Info", "pass")).encode())

        # confirm that we get pdu line start, indicating that connection was successful
        res = self.telnet.expect([bytes("pdu#0>", "utf-8")], self.tmt)
        if(res[0] == -1):
            raise ConnectionError("No PDU prompt")

        return self

    def __exit__(self, type, value, tb):
        """Closes connection"""

        self.telnet.close()
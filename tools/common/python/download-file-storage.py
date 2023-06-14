# TOOL download-file-storage.py: "Download large file from URL" (Use this tool if the size of the file exceeds 200 GB. Download a file from a URL address to the Chipster server. The URL must be visible to Chipster server. If it's not, download the file first to your laptop and upload it then to Chipster.)
# OUTPUT downloaded_file: "Downloaded file"
# PARAMETER url_str: "URL" TYPE UNCHECKED_STRING (URL to download)
# PARAMETER OPTIONAL file_extension: "Add a file extension" TYPE [current: "Keep current", bam: "BAM", fa: "FASTA", fastq: "FASTQ", gtf: "GTF"] DEFAULT current (The output file is named according to the last part of the URL. If it doesn't contain a correct file extension, select it here so that the file type is recognized correctly.)
# PARAMETER OPTIONAL check_certs_str: "Require valid TLS certificate" TYPE [yes: "Yes", no: "No"] DEFAULT yes (Disable if the https server has a self-signed ssl certificate which you trust.)
# RUNTIME python3
# STORAGE 1000
# TOOLS_BIN ""

# add the tools dir to path, because __main__ script cannot use relative imports
import tool_utils
from urllib.parse import urlparse
import ipaddress
import socket
import requests
from http.client import responses
import time


check_certs = check_certs_str == "yes"

if not check_certs:
	print("certificate verification disabled")

try:
	url = urlparse(url_str)

except ValueError:
	raise RuntimeError("CHIPSTER-NOTE: " + "Invalid url: " + url_str)

allowed_protocols = ["http", "https", "ftp"]

if url.scheme not in allowed_protocols:
	raise RuntimeError("CHIPSTER-NOTE: " + "Unsupported protocol: " + url.scheme)

# seems to work even if the netloc is IP already
netloc_ip_str = socket.gethostbyname(url.netloc)

netloc_ip = ipaddress.ip_address(netloc_ip_str)

# check that the ip address isn't a localhost or private address
# this check would be more reliable, if we would do the donwload request with the 
# this ip address, but that would break virtual hosts
if netloc_ip.is_loopback:
	raise RuntimeError("CHIPSTER-NOTE: " +  "Not allowed to connect to localhost: " + url.netloc)

if netloc_ip.is_private:
	raise RuntimeError("CHIPSTER-NOTE: " +  "Not allowed to connect to private address: " + url.netloc)

# parse filename
fragment_removed = url.path.split("#")[0]  # keep to left of first #
query_string_removed = fragment_removed.split("?")[0]
filename = query_string_removed[query_string_removed.rfind("/") + 1:]

if len(filename) == 0:
	filename = "downloaded_file"

if file_extension != "current":
	filename = filename + "." + file_extension

#TDOO allow downloads with unknown size

response = requests.get(url_str, stream=True, verify=check_certs)

if response.status_code >= 400:
	status_name = responses[response.status_code]
	raise RuntimeError("CHIPSTER-NOTE: HTTP error " + str(response.status_code) + " " + status_name + "\n\n" + response.text)

# print(response.headers.keys())

total_size_in_bytes = int(response.headers.get('Content-Length', -1))

if total_size_in_bytes == -1:
	print("size information is not available")

# Doesn't seem to be very sensitive, at least with ping ~2 ms and speed ~300 MB/s.
# Less than 10 kiB slows down the download, more than 100 MiB affects 
# reporting interval and memory consumption.
block_size = 1024 * 1024 # 1 MiB

downloaded_size = 0
t1 = time.perf_counter()
last_reported_size = 0
start_time = t1

# report every second at first and then gradually raise the interval
report_interval = 1
report_interval_multiplier = 1.2

print("downloading...")

with open('downloaded_file', 'wb') as file:
	for data in response.iter_content(block_size):
		
		downloaded_size += len(data)		

		t2 = time.perf_counter()
		
		if t2 - t1 > report_interval:

			report_interval = report_interval * report_interval_multiplier
			speed = (downloaded_size - last_reported_size) / (t2 - t1)
			t1 = t2
			last_reported_size = downloaded_size		

			print(tool_utils.human_readable(downloaded_size) + " / " + tool_utils.human_readable(total_size_in_bytes) + " \t" + tool_utils.human_readable(speed) + "/s")			

		file.write(data)

if total_size_in_bytes != -1 and downloaded_size != total_size_in_bytes:
	raise RuntimeError("CHIPSTER-NOTE: " + "Download failed. Downloaded size " + tool_utils.human_readable(downloaded_size) + " does not match size reported by the server " + tool_utils.human_readable(total_size_in_bytes))

tool_utils.write_output_definitions({"downloaded_file": filename})

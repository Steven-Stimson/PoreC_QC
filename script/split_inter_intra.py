#!/usr/bin/env python3
import sys
import os
import gzip
import threading
from queue import Queue

BUFFER_SIZE = 1_000_000   # 1 MB chunks for faster I/O
QUEUE_SIZE  = 100         # avoid memory blowup

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

# Thread worker writing chunks to file
def writer_worker(out_path, q):
    with gzip.open(out_path, "wt", compresslevel=4) as fout:
        while True:
            chunk = q.get()
            if chunk is None:   # poison pill
                break
            fout.write(chunk)
            q.task_done()
        q.task_done()

if len(sys.argv) != 3:
    sys.exit(f"python3 {sys.argv[0]} <pairs> <outDir>")

pairFile = sys.argv[1]
outDir = sys.argv[2]

baseName = os.path.basename(pairFile).replace(".gz", "")
outFile_intra = os.path.join(outDir, baseName + ".intra.pairs.gz")
outFile_inter = os.path.join(outDir, baseName + ".inter.pairs.gz")

# Queues for multithreaded write
q_intra = Queue(maxsize=QUEUE_SIZE)
q_inter = Queue(maxsize=QUEUE_SIZE)

# Launch writer threads
t_intra = threading.Thread(target=writer_worker, args=(outFile_intra, q_intra), daemon=True)
t_inter = threading.Thread(target=writer_worker, args=(outFile_inter, q_inter), daemon=True)
t_intra.start()
t_inter.start()

buf_intra = []
buf_inter = []
cur_size_intra = 0
cur_size_inter = 0

with open_maybe_gz(pairFile) as f:
    for line in f:
        if line.startswith("#"):
            # Header lines → both outputs
            buf_intra.append(line)
            buf_inter.append(line)
            cur_size_intra += len(line)
            cur_size_inter += len(line)
        else:
            tmp = line.split("\t")
            if tmp[1] == tmp[3]:
                buf_intra.append(line)
                cur_size_intra += len(line)
            else:
                buf_inter.append(line)
                cur_size_inter += len(line)

        # Flush buffers when large
        if cur_size_intra >= BUFFER_SIZE:
            q_intra.put("".join(buf_intra))
            buf_intra = []
            cur_size_intra = 0

        if cur_size_inter >= BUFFER_SIZE:
            q_inter.put("".join(buf_inter))
            buf_inter = []
            cur_size_inter = 0

# Flush any remaining data
if buf_intra:
    q_intra.put("".join(buf_intra))
if buf_inter:
    q_inter.put("".join(buf_inter))

# Send poison pill
q_intra.put(None)
q_inter.put(None)

# Wait for writers
q_intra.join()
q_inter.join()

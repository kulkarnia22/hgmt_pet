import struct
file_path = "../data/debug4.data"

"""with open(file_path, "rb") as f:
    for k in range(100):  # Just read a few entries
        raw = f.read(4)
        if len(raw) < 4:
            print("EOF reached too early")
            break
        i = struct.unpack("<i", raw)[0]
        print(f"int #{k}: {i}")"""

annihilations = []
current = []

with open(file_path, "rb") as f:
    while True:
        header = f.read(4)
        if len(header) < 4:
            print("header length < 4")
            break  # End of file

        i = struct.unpack("<i", header)[0]

        if i == -1:
            print("finished annihilation " + str(len(annihilations)))
            if current:
                annihilations.append(current)
                current = []
            if len(annihilations) == 40:
                print(annihilations[39])
                break
                #don't care about other data for now
        else:
            j_bytes = f.read(4)
            d_bytes = f.read(8)
            if len(j_bytes) < 4 or len(d_bytes) < 8:
                print("j_byte < 4 or d_byte < 8")
                break  # Not enough data left for a full record
            j = struct.unpack("<i", j_bytes)[0]
            impact = struct.unpack("<d", d_bytes)[0]
            current.append((i, j, impact))
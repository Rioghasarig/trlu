import os
matfiles = os.listdir('test_data');
matfiles = [m for m in matfiles if m[-4:] == '.mat']
for i, filename in enumerate(matfiles):
    matname = filename[:-4]
    callgrind_out = f"callgrind/callgrind.{matname}"
    callgrind_options = f"--callgrind-out-file={callgrind_out} --dump-after=clusol"
    valgrind = f"valgrind --tool=callgrind {callgrind_options}"
    ml_script = f"factor('test_data/{filename}'); exit;"
    command = f"""matlab -nojvm -r "{ml_script}" -D"{valgrind}" """
    print(f"{i}: {matname}")
    os.system(command)
    

import subprocess

lrsnash_expected = """dec_digits=1000  lrs_digits=112

2  1/3  2/3  4 
1  2/3  1/3  0  2/3 

2  2/3  1/3  3 
1  0  1/3  2/3  8/3 

2  1  0  3 
1  0  0  1  4 

*Number of equilibria found: 3
*Player 1: vertices=5 bases=5 pivots=8
*Player 2: vertices=3 bases=1 pivots=8"""

lrsnash1_expected = """*lrsnash::lrslib_v.7.2_2022.3.6(64bit,lrslong.h)
2  1/3  2/3  4 
1  2/3  1/3  0  2/3 

2  2/3  1/3  3 
1  0  1/3  2/3  8/3 

2  1  0  3 
1  0  0  1  4 

*Number of equilibria found: 3
*Player 1: vertices=5 bases=5 pivots=8
*Player 2: vertices=3 bases=1 pivots=8

*overflow checking on lrslong arithmetic
*lrsnash::lrslib_v.7.2_2022.3.6(64bit,lrslong.h)"""

lrsnash2_expected = """*lrsnash::lrslib_v.7.2_2022.3.6(128bit,lrslong.h)
2  1/3  2/3  4 
1  2/3  1/3  0  2/3 

2  2/3  1/3  3 
1  0  1/3  2/3  8/3 

2  1  0  3 
1  0  0  1  4 

*Number of equilibria found: 3
*Player 1: vertices=5 bases=5 pivots=8
*Player 2: vertices=3 bases=1 pivots=8

*overflow checking on lrslong arithmetic
*lrsnash::lrslib_v.7.2_2022.3.6(128bit,lrslong.h)"""

lrsnashgmp_expected = """*lrsnash::lrslib_v.7.2_2022.3.6(64bit,lrsgmp.h)_gmp_v.6.3
2  1/3  2/3  4 
1  2/3  1/3  0  2/3 

2  2/3  1/3  3 
1  0  1/3  2/3  8/3 

2  1  0  3 
1  0  0  1  4 

*Number of equilibria found: 3
*Player 1: vertices=5 bases=5 pivots=8
*Player 2: vertices=3 bases=1 pivots=8

*lrsnash::lrslib_v.7.2_2022.3.6(64bit,lrsgmp.h)"""

def call_executable_and_compare(process_name, game_path, expected_output, quiet):
    try:
        # Replace 'executable_name' with the actual executable you want to call
        result = subprocess.run([process_name, game_path], capture_output=True, text=True, check=True)
        result = result.stdout.strip()  # Get stdout and strip any surrounding whitespace
        lines = result.split('\n')[:] if quiet  else result.split('\n')[:-1]
        result = '\n'.join(lines);
        if result == expected_output:
            print("Output matches the expected string.")
        else:
            print("Output does not match the expected string.")
            print(f"Expected: {expected_output}")
            print(f"Actual  : {result}")
    except subprocess.CalledProcessError as e:
        print(f"Error calling executable: {e}")


call_executable_and_compare("./build/lrsnash", "./games/game", lrsnash_expected, True)
call_executable_and_compare("./build/lrsnash1", "./games/game", lrsnash1_expected, False)
call_executable_and_compare("./build/lrsnash2", "./games/game", lrsnash2_expected, False)
call_executable_and_compare("./build/lrsnashgmp", "./games/game", lrsnashgmp_expected, False)

call_executable_and_compare("./lrsnash", "./games/game", lrsnash_expected, True)
call_executable_and_compare("./lrsnash1", "./games/game", lrsnash1_expected, False)
call_executable_and_compare("./lrsnash2", "./games/game", lrsnash2_expected, False)
call_executable_and_compare("./lrsnashgmp", "./games/game", lrsnashgmp_expected, False)

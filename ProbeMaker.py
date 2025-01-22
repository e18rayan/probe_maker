import re

def polya_trim(sequence) -> str:
    '''Trims the poly(A) tail from sequence
        assumes 3'- position of tail'''
    sequence = sequence.upper()
    return sequence.rstrip("A")

def gc_content(probe_test) -> bool:
    '''Calculates the GC content of the input sequence
        returns True if it's in the desired 40-60% range '''
    num_A = 0
    num_T = 0
    num_C = 0
    num_G = 0
    num_N = 0
    for nuc in probe_test:
        match nuc:
            case "A":
                num_A += 1
            case "T":
                num_T += 1
            case "C":
                num_C += 1
            case "G":
                num_G += 1
            case _:
                num_N += 1
    if 0.4 <= (num_C + num_G)/(num_C + num_G + num_A + num_T + num_N) <= 0.6:
        return True
    else:
        return False 

def homopolymer_test(probe_test2) -> bool:
    '''Does a regex search for any nucleotide that is longer than 5 nucleotides
        returns True if it does, otherwise it's okay for next step'''
    if re.search(r'([ATCG])\1{6,}', probe_test2):
        return True
    else:
        return False

def get_probes(sequence, len_probes, overlap) -> list:
    '''Requires 3 inputs, the sequence, the desired length of probes,
        and desired overlap of probes.
        Starts at position 0 and goes all the way until the length of probes
        It restarts by 1 position if the probe fails at GC content or homopolymer.
        Once a probe is added, the loop will subtract the overlap from the final posiiton,
        then continue from that point on to make another probe until the end'''
    probe_sequences = []
    starting_point = 0
    
    while True:
        ending_point = starting_point + int(len_probes)
        if ending_point > len(sequence):
            #breaks if it cannot make another probe of desired length
            break
        else:
            probe_n = sequence[starting_point:ending_point]
            if gc_content(probe_n):
                #checks for GC content between 40%-60%
                if not homopolymer_test(probe_n):
                    #checks for homopolymers, then adds it if it fails
                    probe_sequences.append(probe_n)
                    starting_point += (len_probes - overlap)
                else:
                    starting_point += 1
            else:
                starting_point += 1
        
    return probe_sequences
  
if __name__ == "__main__":
    #No checks if invalid characters or integers
    input_sequence = input(str("What's the sequence?: ")).strip().upper()
    input_length = int(input("How long are the probes?: ").strip())
    input_overlap = int(input("Overlap of probes?: ").strip())
    input_sequence = polya_trim(input_sequence)
    desired_probes = get_probes(input_sequence, input_length, input_overlap)
    
    for pos, seq in enumerate(desired_probes):
        print(f"Probe {pos+1}: {seq}")

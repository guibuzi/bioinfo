import sys
import s01_initialize as init


if __name__ == "__main__":
    data = init.read_fasta("mosaic/us_ha_sampling.fasta")
    seq_list = list(data.values())
    popu = init.main(seq_list, 4, len(seq_list))
    print(popu[0][0])
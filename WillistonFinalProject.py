#1. Read from XML file - similar to that of the optional HW problem
#2. Pull out info for particular spectrum
#3. Given peptide sequence to annotate spectrum with
#4. X-axis of mass or mass-to-charge ratio
#5. Y-axis of abundance
#6. Vertical bars = peaks of spectrum and are fragments of peptide
#7. Y-ions number fragments from C-terminus of peptide
#8. B-ions number fragments from N-terminus of peptide
#9. bi = sum(molecular weight) + 1 == charge state 1
#10. yi = sum(molecular weight of (n+1 - j) + 19)
#11. Molecular Weight of Glycine = 57 DO NOT USE 75
#12. Given peptide, calculate MW of b-ions and y-ions
#13. Look at peaks in spectrum which you read from XML file
#14. Figure out which ones match
#15. DONT EXPECT TO SEE ALL IONS MATCH
#16. DONT EXPECT TO MATCH ALL PEAKS


import sys
from base64 import b64decode
from array import array
import xml.etree.ElementTree as ET
import gzip
import matplotlib.pyplot as plt

#requires 3 command-line arguments
if len(sys.argv) != 4:
    print("Usage: python WillistonFinalProject.py 17mix_test2.mzxml.gz <scan_number> <peptide_seq>")
    sys.exit(1)

#try opening the file, handle errors related to file access
try:
    xml = gzip.open(sys.argv[1])
except FileNotFoundError:
    print("Error: The file", sys.argv[1], "was not found.")
    sys.exit(1)

#validate scan number input
try:
    target_scan = int(sys.argv[2])
except ValueError:
    print("Error: Scan number must be an integer.")
    sys.exit(1)

#validate peptide seq input
peptide_seq = sys.argv[3].upper()
valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
#check if the peptide sequence contains only 20 valid amino acids
invalid_aa_count = 0
invalid_aa_list = []
for aa in peptide_seq:
    if aa not in valid_amino_acids:
        invalid_aa_count += 1
        invalid_aa_list.append(aa)
if invalid_aa_count > 0:
    print("Warning:", invalid_aa_count, "invalid amino acid(s) detected:", ', '.join(invalid_aa_list))
    sys.exit(1)


#dictionary mapping amino acids to its MW in Da
amino_acid_mw = {
    'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04,
    'F': 147.07, 'G': 57.02, 'H': 137.06, 'I': 113.08,
    'K': 128.09, 'L': 113.08, 'M': 131.04, 'N': 114.04,
    'P': 97.05, 'Q': 128.06, 'R': 156.10, 'S': 87.03,
    'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06
}


#extracts peaks for specific scan from XML
def peaks(xml_file, target_scan):
    ns = "{http://sashimi.sourceforge.net/schema/}"
    try:
        context = ET.iterparse(xml_file, events=('start', 'end'))
    except Exception as e:
        print("Error parsing XML:", e)
        return None

    for event, elem in context:
        if event == 'start' and elem.tag == ns + 'scan':
            scan_num = elem.get('num')
            if scan_num == str(target_scan):
                peakselt = elem.find(ns + 'peaks')
                if peakselt is not None:
                    #print("Debug: Found <peaks> for scan", target_scan)
                    #print("Debug: peakselt.text = ",peakselt.text)
                    if peakselt.text is None or not peakselt.text.strip():
                        print("Warning: <peaks> is empty for scan", target_scan)
                        return None
                    try:
                        peaks = array('f', b64decode(peakselt.text))
                        if sys.byteorder != 'big':
                            peaks.byteswap()
                        mzs = peaks[::2]
                        ints = peaks[1::2]
                        return list(zip(mzs, ints))  #return list of (m/z, intensity) pairs
                    except Exception as e:
                        print("Error decoding peaks data for scan number", target_scan,":",e)
                        return None
        elem.clear()
    print("Error: Scan number", target_scan, "not found in the file.")
    return None


#calcs m/z values of b-ions and y-ions for peptide seq
def calc_ions(peptide):
    if not peptide:
        print("Error: Peptide sequence is empty or invalid.")
        return [], []

    b_ions_mz = []
    y_ions_mz = []
    peptide_length = len(peptide)
    b_mass = 0
    y_mass = 0

    #sums MW of amino acids
    for i in range(peptide_length):
        if peptide[i] not in amino_acid_mw:
            print("Error: Invalid amino acid", peptide[i], "in peptide sequence.")
            return [], []
        b_mass += amino_acid_mw[peptide[i]]
        b_ions_mz.append((b_mass + 1) / 1)
        y_mass += amino_acid_mw[peptide[peptide_length - 1 - i]]
        y_ions_mz.append((y_mass + 19) / 1)

    if not b_ions_mz or not y_ions_mz:
        print("Error: No valid b-ions or y-ions calculated for the peptide sequence.")
        
    return b_ions_mz, y_ions_mz


#calcs m/z values of b-ions and y-ions for peptide seq
def plt_ions(mz_values, intensities, b_ions, y_ions, intensity_threshold=0.05):
    if not mz_values or not intensities:
        print("Error: No spectrum data available to plot.")
        return

    max_intensity = max(intensities) if intensities else 1
    relative_abundance = [(i / max_intensity) * 100 for i in intensities]
    #print("Max Intensity:", max_intensity)
    #print("Relative Abundance Range:",min(relative_abundance),max(relative_abundance))

    matched_b_ions = []
    matched_y_ions = []
    unmatched_peaks = []


    #match m/z values from spectrum to b-ions and y-ions
    for mz, rel_ab in zip(mz_values, relative_abundance):
        if rel_ab < intensity_threshold * 100:
            continue

        #includes a tolerance
        matched_b = any(abs(mz-b) < 0.05 for b in b_ions)
        matched_y = any(abs(mz-y) < 0.05 for y in y_ions)

        if matched_b:
            matched_b_ions.append((mz, rel_ab))
        elif matched_y:
            matched_y_ions.append((mz, rel_ab))
        else:
            unmatched_peaks.append((mz, rel_ab))

    if matched_b_ions or matched_y_ions:
        # Print matched b-ions and y-ions with m/z, index, and relative abundance
        if matched_b_ions:
            print("Matched b-ions:")
            for i, (mz, rel_ab) in enumerate(matched_b_ions, start=1):
                print("b"+str(i),": m/z =", round(mz,2), "Relative Abundance =", round(rel_ab,2),"%")

        if matched_y_ions:
            print("Matched y-ions:")
            for i, (mz, rel_ab) in enumerate(matched_y_ions, start=1):
                print("y"+str(i),": m/z =", round(mz,2), "Relative Abundance =", round(rel_ab,2),"%")
    else:
        print("Warning: No b-ions or y_ions found in spectrum for peptide", peptide_seq)


    #plot unmatched peaks
    plt.figure(figsize=(12, 6))
    if mz_values:
        plt.stem(mz_values, relative_abundance, linefmt='k-', markerfmt=' ', basefmt=' ', label='Unmatched Peaks')

    #plot matched b-ions and y-ions
    if matched_b_ions:
        mzs, abundance = zip(*matched_b_ions)
        plt.stem(mzs, abundance, linefmt='blue', markerfmt=' ', basefmt=' ', label='b-ions')

    if matched_y_ions:
        mzs, abundance = zip(*matched_y_ions)
        plt.stem(mzs, abundance, linefmt='red', markerfmt=' ', basefmt=' ', label='y-ions')

    #annotate b-ions and y-ions
    if matched_b_ions:
        for i, (mz, rel_ab) in enumerate(matched_b_ions, start=1):
            plt.text(mz, rel_ab + 3, 'b' + str(i), color='blue', ha='center', fontsize=10)

    if matched_y_ions:
        for i, (mz, rel_ab) in enumerate(matched_y_ions, start=1):
            plt.text(mz, rel_ab + 3, 'y' + str(i), color='red', ha='center', fontsize=10)

    #set axis limits
    x_min = min(mz_values) - 10
    x_max = max(mz_values) + 10
    y_max = max(relative_abundance) + 5

    plt.xlabel("Mass to Charge Ratio (m/z)", fontsize=12)
    plt.ylabel("Relative Abundance (%)", fontsize=12)
    plt.title("Spectrum with Annotated b-ions and y-ions", fontsize=14)
    plt.ylim(0, y_max)
    plt.xlim(x_min, x_max)
    plt.legend(loc='upper right', fontsize=10)
    plt.show()


#extract peaks data for specified scan number
peaks_data = peaks(xml, target_scan)
if peaks_data:
    mz_values, intensities = zip(*peaks_data)
    b_ions, y_ions = calc_ions(peptide_seq)
    #print("B-Ions:",b_ions)
    #print("Y-Ions:",y_ions)
    if b_ions and y_ions:
        plt_ions(mz_values, intensities, b_ions, y_ions)
    else:
        print("Error: Could not calculate b_ions or y_ions for the given peptide")
        sys.exit(1)
else:
    print("No peaks found in given scan")
    sys.exit(1)

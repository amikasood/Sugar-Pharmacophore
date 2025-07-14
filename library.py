import collections

class MonosaccharideAnalyzer:
    def __init__(self):
        pass

    def _parse_ehf(self, ehf_string):
        """
        Parses an EHF (Extended Haworth Format) string into a structured representation.
        Handles ring flipping internally for consistency.
        """
        parts = ehf_string.strip().split('_')
        if len(parts) != 6:
            raise ValueError("EHF string must represent exactly 6 ring atoms.")

        monosaccharide = {'bb': [None] * 12, 'eq': [None] * 12, 'ax': [None] * 12}

        for i, part in enumerate(parts):
            sub_parts = part.split('-')
            monosaccharide['bb'][i] = sub_parts[0]

            if len(sub_parts) > 1:
                # Determine if it's axial (lowercase) or equatorial (uppercase)
                if sub_parts[1][0].islower():
                    monosaccharide['ax'][i] = sub_parts[1]
                else:
                    monosaccharide['eq'][i] = sub_parts[1]

            if len(sub_parts) > 2:
                if sub_parts[2][0].islower():
                    if monosaccharide['ax'][i] is not None:
                        raise ValueError(f"More than one axial group at position {i+1}")
                    monosaccharide['ax'][i] = sub_parts[2]
                else:
                    if monosaccharide['eq'][i] is not None:
                        raise ValueError(f"More than one equatorial group at position {i+1}")
                    monosaccharide['eq'][i] = sub_parts[2]

            # Fill in empty positions with a space
            if monosaccharide['eq'][i] is None:
                monosaccharide['eq'][i] = " "
            if monosaccharide['ax'][i] is None:
                monosaccharide['ax'][i] = " "
            
            # Validate ring atom
            if monosaccharide['bb'][i][0].lower() not in ['c', 'o']:
                raise ValueError(f"Invalid ring atom '{monosaccharide['bb'][i][0]}' at position {i+1}")


        # Create flipped version for carbons and oxygens only
        flipped_monosaccharide = {'bb': [None] * 12, 'eq': [None] * 12, 'ax': [None] * 12}
        
        # Ring O in front becomes Ring O in back.
        # C1 becomes C5, C2 becomes C4, C3 remains C3, C4 becomes C2, and C5 becomes C1.
        # The logic is (original_index + offset) % 6
        # This means 0 -> 0 (O), 1 -> 5 (C5), 2 -> 4 (C4), 3 -> 3 (C3), 4 -> 2 (C2), 5 -> 1 (C1)
        # Essentially, the order of atoms is reversed around the O-C3 axis.
        # Original: bb[0], bb[1], bb[2], bb[3], bb[4], bb[5]
        # Flipped:  bb[0], bb[5], bb[4], bb[3], bb[2], bb[1]
        
        # Position 0 (O) remains in place, but its case is swapped if it's C or O
        flipped_monosaccharide['bb'][0] = monosaccharide['bb'][0].swapcase() if monosaccharide['bb'][0].lower() in ['c','o'] else monosaccharide['bb'][0]
        flipped_monosaccharide['ax'][0] = monosaccharide['ax'][0] # Axial/equatorial at O is usually empty or irrelevant
        flipped_monosaccharide['eq'][0] = monosaccharide['eq'][0]
        
        # Remaining positions are flipped symmetrically
        for j in range(1, 6):
            flipped_monosaccharide['ax'][6 - j] = monosaccharide['ax'][j]
            flipped_monosaccharide['eq'][6 - j] = monosaccharide['eq'][j]
            flipped_monosaccharide['bb'][6 - j] = monosaccharide['bb'][j].swapcase() if monosaccharide['bb'][j].lower() in ['c','o'] else monosaccharide['bb'][j]


        # Copy contents for positions 6-11 (for cyclic comparison)
        for j in range(6, 12):
            monosaccharide['bb'][j] = monosaccharide['bb'][j - 6]
            monosaccharide['ax'][j] = monosaccharide['ax'][j - 6]
            monosaccharide['eq'][j] = monosaccharide['eq'][j - 6]
            flipped_monosaccharide['bb'][j] = flipped_monosaccharide['bb'][j - 6]
            flipped_monosaccharide['ax'][j] = flipped_monosaccharide['ax'][j - 6]
            flipped_monosaccharide['eq'][j] = flipped_monosaccharide['eq'][j - 6]
        return monosaccharide, flipped_monosaccharide


    def score_structural_similarity(self, ehf_string1, ehf_string2):
        """
        Compares two monosaccharide EHF strings and returns their structural similarity score.
        
        Args:
            ehf_string1 (str): EHF representation of the first monosaccharide.
            ehf_string2 (str): EHF representation of the second monosaccharide.

        Returns:
            float: The best (highest) similarity score among all possible orientations.
            
        Raises:
            ValueError: If the EHF string format is incorrect.
        """
        mono1_normal, mono1_flipped = self._parse_ehf(ehf_string1)
        mono2_normal, mono2_flipped = self._parse_ehf(ehf_string2)

        def calculate_score(m1, m2):
            best_score = -100
            
            for k in range(6):  # Starting position for alignment in m1
                current_score = 0
                
                for l in range(6):  # Position in m2
                    # Check ring atoms
                    if m1['bb'][k + l] == m2['bb'][l]:
                        current_score += 2.0  # Same plane and same atom
                        if m1['bb'][k + l].lower() == 'c':
                            if m1['ax'][k + l] == m2['ax'][l]:
                                current_score += 1
                            if m1['eq'][k + l] == m2['eq'][l]:
                                current_score += 1
                        elif m1['bb'][k+l].lower() == 'o':
                            current_score += 2.0
                    else:
                        # Ring atoms don't match, check if they are in the same plane
                        if (m1['bb'][k + l].lower() == 'o' and m2['bb'][l].lower() == 'c') or \
                           (m1['bb'][k + l].lower() == 'c' and m2['bb'][l].lower() == 'o'):
                            current_score += 1.0
                            if m1['ax'][k + l] == m2['ax'][l]:
                                current_score += 1
                            if m1['eq'][k + l] == m2['eq'][l]:
                                current_score += 1
                
                best_score = max(best_score, current_score)
            return best_score
        
        scores = []
        scores.append(calculate_score(mono1_normal, mono2_normal))
        scores.append(calculate_score(mono1_normal, mono2_flipped))
        scores.append(calculate_score(mono1_flipped, mono2_normal))
        scores.append(calculate_score(mono1_flipped, mono2_flipped))
        
        return max(scores)


    def predict_pharmacophore(self, ehf_strings):
        """
        Predicts a common pharmacophore from a list of monosaccharide EHF strings.
        
        Args:
            ehf_strings (list): A list of EHF representation strings of monosaccharides.
                It is recommended that the structures are listed in descending order of affinity.

        Returns:
            dict: A dictionary representing the predicted pharmacophore, with 'bb', 'ax', and 'eq' lists,
                  and the relative scores for each position.
            
        Raises:
            ValueError: If the input list is empty or contains malformed EHF strings.
        """
        if not ehf_strings:
            raise ValueError("Input list of EHF strings cannot be empty.")
        if len(ehf_strings) == 1:
            print("\n\tThere is only one structure in the input file. Pharmacophore prediction requires multiple structures.\n")
            return None

        # Parse all monosaccharides and their flipped versions
        monosaccs_raw = []
        for s in ehf_strings:
            normal, flipped = self._parse_ehf(s)
            monosaccs_raw.append(normal)
            monosaccs_raw.append(flipped) # Append flipped structures as well

        nstr = len(ehf_strings)
        aligned_monosaccs = [{} for _ in range(nstr)]

        # Align all structures to the first structure (monosaccs_raw[0] - normal orientation of first input)
        # The logic for alignment is similar to score_structural_similarity but it keeps the best aligned structure
        
        # Initialize aligned_monosaccs[0] with the first structure
        aligned_monosaccs[0]['bb'] = monosaccs_raw[0]['bb'][0:6]
        aligned_monosaccs[0]['ax'] = monosaccs_raw[0]['ax'][0:6]
        aligned_monosaccs[0]['eq'] = monosaccs_raw[0]['eq'][0:6]

        for i in range(1, nstr):
            best_score = -100
            best_aligned_structure = {'bb': [None]*6, 'ax': [None]*6, 'eq': [None]*6} # Initialize here

            # Check normal orientation of current monosaccharide (monosaccs_raw[i])
            for k in range(6): # offset for alignment
                current_score = 0
                for l in range(6):
                    current_bb = monosaccs_raw[i]['bb'][k+l]
                    current_ax = monosaccs_raw[i]['ax'][k+l]
                    current_eq = monosaccs_raw[i]['eq'][k+l]
                    
                    target_bb = monosaccs_raw[0]['bb'][l]
                    target_ax = monosaccs_raw[0]['ax'][l]
                    target_eq = monosaccs_raw[0]['eq'][l]

                    if current_bb == target_bb:
                        current_score += 2.0
                        if current_bb.lower() == 'c':
                            if current_ax == target_ax:
                                current_score += 1
                            if current_eq == target_eq:
                                current_score += 1
                        elif current_bb.lower() == 'o':
                            current_score += 2.0
                    else:
                        if (current_bb.lower() == 'o' and target_bb.lower() == 'c') or \
                           (current_bb.lower() == 'c' and target_bb.lower() == 'o'):
                            current_score += 1.0
                            if current_ax == target_ax:
                                current_score += 1
                            if current_eq == target_eq:
                                current_score += 1
                
                if current_score > best_score:
                    best_score = current_score
                    best_aligned_structure['bb'] = monosaccs_raw[i]['bb'][k:k+6]
                    best_aligned_structure['ax'] = monosaccs_raw[i]['ax'][k:k+6]
                    best_aligned_structure['eq'] = monosaccs_raw[i]['eq'][k:k+6]
            
            # Check flipped orientation of current monosaccharide (monosaccs_raw[i + nstr])
            for k in range(6): # offset for alignment
                current_score = 0
                flipped_index = i + nstr 
                for l in range(6):
                    current_bb = monosaccs_raw[flipped_index]['bb'][k+l]
                    current_ax = monosaccs_raw[flipped_index]['ax'][k+l]
                    current_eq = monosaccs_raw[flipped_index]['eq'][k+l]
                    
                    target_bb = monosaccs_raw[0]['bb'][l]
                    target_ax = monosaccs_raw[0]['ax'][l]
                    target_eq = monosaccs_raw[0]['eq'][l]

                    if current_bb == target_bb:
                        current_score += 2.0
                        if current_bb.lower() == 'c':
                            if current_ax == target_ax:
                                current_score += 1
                            if current_eq == target_eq:
                                current_score += 1
                        elif current_bb.lower() == 'o':
                            current_score += 2.0
                    else:
                        if (current_bb.lower() == 'o' and target_bb.lower() == 'c') or \
                           (current_bb.lower() == 'c' and target_bb.lower() == 'o'):
                            current_score += 1.0
                            if current_ax == target_ax:
                                current_score += 1
                            if current_eq == target_eq:
                                current_score += 1

                if current_score > best_score:
                    best_score = current_score
                    best_aligned_structure['bb'] = monosaccs_raw[flipped_index]['bb'][k:k+6]
                    best_aligned_structure['ax'] = monosaccs_raw[flipped_index]['ax'][k:k+6]
                    best_aligned_structure['eq'] = monosaccs_raw[flipped_index]['eq'][k:k+6]
            
            aligned_monosaccs[i] = best_aligned_structure
            print(f"alugned str: i, {aligned_monosaccs[i]}\n")

        # Calculate unique scores and determine pharmacophore
        unique_scores = {
            'bb': [collections.defaultdict(int) for _ in range(6)],
            'ax': [collections.defaultdict(int) for _ in range(6)],
            'eq': [collections.defaultdict(int) for _ in range(6)]
        }

        for i in range(nstr):
            for l in range(6):
                unique_scores['bb'][l][aligned_monosaccs[i]['bb'][l]] += 1
                unique_scores['ax'][l][aligned_monosaccs[i]['ax'][l]] += 1
                unique_scores['eq'][l][aligned_monosaccs[i]['eq'][l]] += 1

        pharmacophore = {
            'bb': [None] * 6,
            'ax': [None] * 6,
            'eq': [None] * 6,
            'scores': {'bb': [0] * 6, 'ax': [0] * 6, 'eq': [0] * 6}
        }

        for l in range(6):
            # Find the most frequent backbone atom
            if unique_scores['bb'][l]:
                pharmacophore['bb'][l] = max(unique_scores['bb'][l], key=unique_scores['bb'][l].get)
                pharmacophore['scores']['bb'][l] = unique_scores['bb'][l][pharmacophore['bb'][l]] / nstr
            else:
                pharmacophore['bb'][l] = " " # No common backbone element found

            # Find the most frequent axial group
            if unique_scores['ax'][l]:
                pharmacophore['ax'][l] = max(unique_scores['ax'][l], key=unique_scores['ax'][l].get)
                pharmacophore['scores']['ax'][l] = unique_scores['ax'][l][pharmacophore['ax'][l]] / nstr
            else:
                pharmacophore['ax'][l] = " " # No common axial group found

            # Find the most frequent equatorial group
            if unique_scores['eq'][l]:
                pharmacophore['eq'][l] = max(unique_scores['eq'][l], key=unique_scores['eq'][l].get)
                pharmacophore['scores']['eq'][l] = unique_scores['eq'][l][pharmacophore['eq'][l]] / nstr
            else:
                pharmacophore['eq'][l] = " " # No common equatorial group found
        
        return pharmacophore


    def find_binders(self, ehf_strings, pharmacophore_string):
        """
        Identifies which monosaccharides from a list contain a given pharmacophore.
        This function is based on the logic in 'binders.c'.
        
        Args:
            ehf_strings (list): A list of EHF representation strings of monosaccharides.
            pharmacophore_string (str): EHF representation of the pharmacophore to search for.

        Returns:
            list: A list of original EHF strings that contain the pharmacophore.
            
        Raises:
            ValueError: If the input list is empty or contains malformed EHF strings or pharmacophore.
        """
        if not ehf_strings:
            raise ValueError("Input list of EHF strings cannot be empty.")
        
        # Parse the pharmacophore string
        # The pharmacophore parsing is slightly different as it doesn't flip itself
        pharm_parts = pharmacophore_string.strip().split('_')
        pharm_size = len(pharm_parts)
        if pharm_size > 6:
            raise ValueError("Pharmacophore string must represent at most 6 ring atoms.")

        pharmacophore = {'bb': [None] * pharm_size, 'eq': [None] * pharm_size, 'ax': [None] * pharm_size}

        for i, part in enumerate(pharm_parts):
            sub_parts = part.split('-')
            pharmacophore['bb'][i] = sub_parts[0]

            if len(sub_parts) > 1:
                if sub_parts[1][0].islower():
                    pharmacophore['ax'][i] = sub_parts[1]
                else:
                    pharmacophore['eq'][i] = sub_parts[1]

            if len(sub_parts) > 2:
                if sub_parts[2][0].islower():
                    if pharmacophore['ax'][i] is not None:
                        raise ValueError(f"More than one axial group in pharmacophore at position {i+1}")
                    pharmacophore['ax'][i] = sub_parts[2]
                else:
                    if pharmacophore['eq'][i] is not None:
                        raise ValueError(f"More than one equatorial group in pharmacophore at position {i+1}")
                    pharmacophore['eq'][i] = sub_parts[2]

            if pharmacophore['eq'][i] is None:
                pharmacophore['eq'][i] = " "
            if pharmacophore['ax'][i] is None:
                pharmacophore['ax'][i] = " "
            
            if pharmacophore['bb'][i][0].lower() not in ['c', 'o', ' ']:
                raise ValueError(f"Invalid ring atom '{pharmacophore['bb'][i][0]}' in pharmacophore at position {i+1}")

        # Parse all monosaccharides and their flipped versions
        monosaccs_processed = []
        original_strings = []
        for s in ehf_strings:
            normal, flipped = self._parse_ehf(s)
            monosaccs_processed.append(normal)
            monosaccs_processed.append(flipped)
            original_strings.append(s)

        nstr = len(ehf_strings)
        matching_structures = []

        for i in range(nstr):
            # Check normal orientation of the monosaccharide
            match_found = False
            for k in range(6):  # Starting position for pharmacophore in monosaccharide
                current_match = True
                for l in range(pharm_size): # Position in pharmacophore
                    
                    if pharmacophore['bb'][l] != " ": # Only compare if pharmacophore has a defined backbone atom
                        if monosaccs_processed[i]['bb'][k+l] != pharmacophore['bb'][l]:
                            current_match = False
                            break
                        
                        if pharmacophore['ax'][l] != " ":
                            if monosaccs_processed[i]['ax'][k+l] != pharmacophore['ax'][l]:
                                current_match = False
                                break
                        
                        if pharmacophore['eq'][l] != " ":
                            if monosaccs_processed[i]['eq'][k+l] != pharmacophore['eq'][l]:
                                current_match = False
                                break
                
                if current_match:
                    match_found = True
                    break
            
            if match_found:
                matching_structures.append(original_strings[i])
                continue # Move to the next monosaccharide

            # If not found in normal orientation, check flipped orientation
            for k in range(6):
                current_match = True
                flipped_index = i + nstr
                for l in range(pharm_size):
                    
                    if pharmacophore['bb'][l] != " ":
                        if monosaccs_processed[flipped_index]['bb'][k+l] != pharmacophore['bb'][l]:
                            current_match = False
                            break
                        
                        if pharmacophore['ax'][l] != " ":
                            if monosaccs_processed[flipped_index]['ax'][k+l] != pharmacophore['ax'][l]:
                                current_match = False
                                break
                        
                        if pharmacophore['eq'][l] != " ":
                            if monosaccs_processed[flipped_index]['eq'][k+l] != pharmacophore['eq'][l]:
                                current_match = False
                                break
                
                if current_match:
                    match_found = True
                    break
            
            if match_found:
                matching_structures.append(original_strings[i])
        
        return matching_structures
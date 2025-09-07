#!/usr/bin/bash

# Global variables are always UPPERCASE.
# Local are used with local keyword and lowercase.
# If some function requires too much arguments,
# try using global variables directly, however, this is harder to
# read and debug.

function ScriptInfo() {
  DATE="2025"
  VERSION="0.0.1"
  GH_URL="https://github.com/tcaceresm/AmberMCPBHelper"

  cat <<EOF
###################################################
 Welcome to prepare4mcpb version ${VERSION} ${DATE}   
  Author: Tomás Cáceres <caceres.tomas@uc.cl>
  Laboratory of Computational simulation & drug design        
  GitHub <${GH_URL}>                             
  Powered by high fat food and procrastination   
###################################################
EOF
}

Help() {
  ScriptInfo
  echo -e "\nUsage: bash prepare4mcpb.sh OPTIONS\n"
  echo "This script parse files for MCPB.py."

  echo "Required options:"
  echo " -d, --work_dir         <path>       Working directory. Inside this directory, a folder named setupMD will be created which contains all necessary files."
  echo " -i, --input_file       <str>        XYZ file."
  echo " --metal_name           <str>        Metal name."
  echo " --metal_charge         <int>        Metal charge (oxidation state)."
  echo " --lig_charge           <int>        Ligand charge (organic fragment)."
  echo " --lig_mult             <int>        Ligand multiplicity."
  echo "Optional:"
  echo " --lig_ff               <str>        (default="gaff2"). Ligand forcefield."
  echo " --lig_charge_method    <str>        (default="abcg2"). Ligand charge method to calculate atom partial charges."
  echo " --lig_resname          <str>        (default="LIG"). Name to rename the residue name of the ligand."
  echo " --ff                   <str>        (default="ff19SB"). Protein forcefield."
  echo " -h, --help                          Show this help."

}

# Check arguments
if [[ "$#" == 0 ]]; then
  echo "Error: No options provided."
  echo "Use --help option to check available options."
  exit 1
fi


# Default values

LIG_RESNAME="LIG"
LIG_FF="gaff2"
LIG_CHARGE_METHOD="abcg2"
PROT_FF="ff19SB"

# CLI option parser
while [[ $# -gt 0 ]]; do
  case "$1" in
  '-d' | '--work_dir'        ) shift ; WDPATH=$1 ;;
  '--input_file' | '-i'      ) shift ; INPUT_FILE=$1 ;;
  '--metal_name'             ) shift ; METAL=$1 ;;
  '--metal_charge'           ) shift ; METAL_CHARGE=$1 ;;
  '--lig_charge'             ) shift ; LIG_CHARGE=$1 ;;
  '--lig_mult'               ) shift ; LIG_MULT=$1 ;;
  '--lig_ff'                 ) shift ; LIG_FF=$1 ;;
  '--lig_charge_method'      ) shift ; LIG_CHARGE_METHOD=$1 ;;
  '--lig_resname'            ) shift ; LIG_RESNAME=$1 ;;
  '--ff'                     ) shift ; PROT_FF=$1 ;;
  '--help' | '-h'            ) Help ; exit 0 ;;
  *                          ) echo "Unrecognized command line option: $1" >> /dev/stderr ; exit 1 ;;
  esac
  shift
done


#### Functions ####

function CheckProgram() {
  for COMMAND in "$@"; do
    # Check if command is available
    if ! command -v ${1} >/dev/null 2>&1; then
      echo "Error: ${1} not available, exiting."
      exit 1
    fi
  done
}

function CheckFiles() {
  # Check existence of files
  for ARG in "$@"; do
    if [[ ! -f ${ARG} ]]; then
      echo "Error: ${ARG} file doesn't exist."
      exit 1
    fi
  done
}

function CheckVariable() {
  # Check if variable is not empty
  for ARG in "$@"; do
    if [[ -z ${ARG} ]]; then
      echo "Error: Required option not provided."
      echo "Use --help option to check available options."
      exit 1
  fi
  done
}

CheckVariable "${WDPATH}" "${INPUT_FILE}" "${METAL}" \
              "${METAL_CHARGE}" "${LIG_RESNAME}" "${LIG_CHARGE}" \
              "${LIG_MULT}" "${LIG_FF}"


function XYZ_2_PDB() {
  # XYZ to PDB format using OpenBabel
  local input_file=$1
  local name=$(basename "${input_file}" .xyz)

  CheckProgram obabel
  CheckFiles ${input_file}

  echo "Converting XYZ file to PDB format"

  obabel -i xyz "${input_file}" -o pdb -O ${name}.pdb || { echo "Error converting file." ; exit 1; }
  
  echo -e "Done.\n"

}

function CleanPDB() {
  # Just keep ATOM|HETATM records in PDB file
  local pdb_file=$1

  CheckFiles "${pdb_file}"

  echo "Cleaning PDB file. Keeping only ATOM or HETATM lines."
  awk '/^(ATOM|HETATM)/' ${pdb_file} > ${pdb_file}.tmp || { echo "Error cleaning PDB." ; exit 1; }
  mv ${pdb_file}.tmp ${pdb_file} 
  echo -e "Done.\n"
}

#function AssignUniqueAtomsName() {
  # Each atom inside a residue must be unique labeled, i.e, atom names must be uniques inside a residue.
  # Input is PDB file.
  # local input_file=$1

  # CheckFiles "${input_file}"

  # input_filename=$(basename ${input_file} .pdb)

  # echo "Assigning unique atoms name..."
  # awk '
  #   BEGIN {
  #   metals["FE"]
  #   metals["RE"]
  #   metals["RU"]
  #   metals["ZN"]
  # } {
  #     elem = toupper(substr($0,13,3))                                 # atom name (col 13-16) https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
  #     gsub(/ /,"",elem)
  #     if ( !(elem in metals) ) {
  #       count[elem]++                                               # associative arrays awk
  #       newname = elem count[elem]                                  # ej: C1, H1, FE1
  #       $0 = substr($0,1,12) sprintf("%-4s", newname) substr($0,17) # substr(string, start, length)
  #     } else {
  #       $0 = substr($0,1,12) sprintf("%-4s", elem) substr($0,17) # substr(string, start, length)  
  #     } 
  #     print
  # }
  # ' "${input_file}" > "${input_filename}_renamed.pdb" || { echo "Error during AssignUniqueAtomsName()." ; exit 1; }
  # echo "Done"

#}

function AssignUniqueAtomsName() {
  # Each atom inside a residue must be unique labeled, i.e, atom names must be uniques inside a residue.
  # Input is PDB file.
  local input_file=$1

  CheckFiles "${input_file}"

  input_filename=$(basename ${input_file} .pdb)

  echo "Assigning unique atoms name..."

  awk -v metal="${METAL}" '
    {
      elem = toupper(substr($0,13,3))     # atom name (col 13-16)
      gsub(/ /,"",elem)                   # eliminar espacios

      if (elem != metal) {
        count[elem]++                      # contador para nombres únicos
        newname = elem count[elem]         # ej: C1, H1, FE1
        $0 = substr($0,1,12) sprintf("%-4s", newname) substr($0,17)
      } else {
        $0 = substr($0,1,12) sprintf("%-4s", elem) substr($0,17)
      }
      print
    }
  
  ' "${input_file}" > "${input_filename}_renamed.pdb" || { echo "Error during AssignUniqueAtomsName()." ; exit 1; }
  
  echo -e "Done.\n"

}

function ExtractMetal() {
  # Extract metal from PDB
  local input_file=$1

  echo "Extracting metal to a new pdb."

  CheckFiles ${input_file}

  # To automate metal identificacion
  # METAL=$(awk '
  #   BEGIN {
  #   metals["FE"]
  #   metals["RE"]
  #   metals["RU"]
  #   metals["ZN"]
  #   }
  #   $3 in metals {
  #     print $3; exit
  #   }
  #   ' "$input_file")

  # Create metal output filename
  local output_filename="${METAL}.pdb"

  # Extract metal from input_file PDB
  awk -v metal="${METAL}" '
    $3 == metal {
      print $0 > metal ".pdb"
    }
  ' "${input_file}" || { echo "Error during extracting metal from pdb."; exit 1; }

  # Global variable. Metal resnumber from initial complex.
  METAL_RESNUMBER=$(awk '{print $5}' ${output_filename})
  METAL_ATOM_ID=$(awk '{print $2}' ${output_filename})

  echo "Renaming resname of created metal pdb."
  local metal_resname=$(awk '{print $4}' ${output_filename})

  sed -i "s/${metal_resname}/${METAL} /g" ${output_filename}
  echo "Created file: ${output_filename}"

  echo "Renumbering metal residue number --> hardcoded to be 1."
  awk '
    {
    $0 = substr($0,1,22) sprintf("%4d", 1) substr($0,27)
    print
    }' ${METAL}.pdb > ${METAL}.pdb.tmp \
  && mv ${METAL}.pdb.tmp ${METAL}.pdb || { echo "Error renumbering metal fragment."; exit 1; }
  echo "Done renumbering."

  echo -e "Done.\n"

}

function ExtractOrganicFragment() {
  local input_file=$1

  CheckFiles ${input_file}

  echo "Extracting organic fragment from ${input_file}"

  awk -v metal=${METAL} -v lig_resname=${LIG_RESNAME} '
  BEGIN {
    output_file = ""
  }
  $3 != metal {
    if (output_file == "") {
      resname = lig_resname
      output_file = resname ".pdb"
    }
    $0 = substr($0,1,17) sprintf("%3s", resname) substr($0,21) # substr(string, start, length)

    print $0 > output_file
  }
  ' "${input_file}" || { echo "Error extracting organic fragment from PDB"; exit 1; }

  echo "Created file: ${LIG_RESNAME}.pdb"
  
  echo "Renumbering ligand residue --> hardcoded to be 2."
  awk '
    {
    $0 = substr($0,1,22) sprintf("%4d", 2) substr($0,27)
    print
    }' ${LIG_RESNAME}.pdb > ${LIG_RESNAME}.pdb.tmp \
  && mv ${LIG_RESNAME}.pdb.tmp ${LIG_RESNAME}.pdb || { echo "Error renumbering organic fragment."; exit 1; }
  echo "Done renumbering."
  
  echo -e "Done.\n"
}

function PrepareMetalMol2() {
  # Get mol2 from metal PDB
  local input_file=$1
  local charge=$2

  CheckFiles ${input_file}
  CheckProgram metalpdb2mol2.py

  echo "Obtaining MOL2 from metal PDB."
  metalpdb2mol2.py -i ${METAL}.pdb -o ${METAL}.mol2 -c ${charge} || { echo "Error obtaining MOL2 from metal PDB" ; exit 1; }  
  echo -e "Done.\n"

}

function PrepareOrganicFragment() {
  # Obtain frcmod file for organic fragment
  local lig_pdb=$1
  local lig_name=$(basename ${lig_pdb} .pdb)

  CheckFiles "${lig_pdb}"
  CheckProgram "parmchk2" "antechamber"

  echo "Running antechamber..."
  antechamber -i ${lig_pdb} -fi pdb -o ${lig_name}.mol2 -fo mol2 -at ${LIG_FF} -c ${LIG_CHARGE_METHOD} -nc ${LIG_CHARGE} -pf y || \
  { echo "Error running antechamber. Exiting"; exit 1; }
  echo -e "Done running antechamber.\n"

  echo "Running parmchk2."
  parmchk2 -i ${lig_name}.mol2 -o ${lig_name}.frcmod -f mol2 || { echo "Error runing parmchk2. Exiting." ; exit 1; }
  echo -e "Done running parmchk2.\n"

}

function LigatingAtoms() {
  # Obtain which atoms coordinates to metal
  # Input is prepared pdb
  local input_file=$1
  local metal_atom_id=$2

  local cutoff=2.5       # distancia de corte en Å
  echo -e "Obtaining ligating atoms..."
  
  if [[ -f "LIGATING_ATOMS.txt" ]]; then
    rm "LIGATING_ATOMS.txt"
  fi

  if [[ -f "PAIRS.txt" ]]; then
    rm "PAIRS.txt"
  fi

  awk -v ref="${metal_atom_id}" -v cutoff="${cutoff}" '
  /^(ATOM|HETATM)/ {
    atom_id = substr($0,7,5)+0
    x = substr($0,31,8)+0
    y = substr($0,39,8)+0
    z = substr($0,47,8)+0

    if (atom_id == ref) 
    {
      xref = x; yref = y; zref = z
    }
    atoms[atom_id,"x"] = x # associative array, like a key-value dictionary
    atoms[atom_id,"y"] = y
    atoms[atom_id,"z"] = z
    lines[atom_id] = $0
  }
  END {

    for (id in lines) # loop through all "keys" which are atoms id in this case
      { 
      if (id == ref) continue # avoid calculating the distance using the metal itself
      dx = atoms[id,"x"] - xref
      dy = atoms[id,"y"] - yref
      dz = atoms[id,"z"] - zref
      d = sqrt(dx*dx + dy*dy + dz*dz) # euclidean distance
      if (d <= cutoff)
        {
        # This format is for MCPB.in file
        printf id "-" ref " " >> "PAIRS.txt"

        # For debugging purposes
        printf("Atom %2d is at %.3f Å -> %s\n", id, d, lines[id]) >> "LIGATING_ATOMS.txt"
        }
      }
  }' ${input_file}

  echo -e "Done.\n"

}

function PrepareMCPBFile() {

  echo "Preparing MCPB input file..."

  if [[ ${LIG_FF} == "gaff2" ]]; then
    local lig_ff=2
  else
    local lig_ff=1
  fi

  CheckVariable "${lig_ff}" "${INPUT_FILENAME}" "${METAL}" \
                "${METAL_ATOM_ID}" "${LIG_RESNAME}" "${PAIRS}" "${PROT_FF}"
  CheckFiles ${INPUT_FILENAME}_prep.pdb ${METAL}.pdb ${LIG_RESNAME}.mol2 ${LIG_RESNAME}.frcmod

  cat > ${LIG_RESNAME}_MCPB.in <<EOF
original_pdb ${INPUT_FILENAME}_prep.pdb
group_name ${INPUT_FILENAME}_prep
cut_off 2.7
ion_ids ${METAL_ATOM_ID}
ion_mol2files ${METAL}.mol2
naa_mol2files ${LIG_RESNAME}.mol2
frcmod_files ${LIG_RESNAME}.frcmod
add_bonded_pairs ${PAIRS}
gaff ${lig_ff}
force_field ${PROT_FF}
software_version g16
smmodel_chg ${LIG_CHARGE}
smmodel_spin ${LIG_MULT}
lgmodel_chg ${LIG_CHARGE}
lgmodel_spin ${LIG_MULT}
EOF
  echo -e "Done.\n"
}


#### End functions ####


#### Main ####

INPUT_FILENAME=$(basename ${INPUT_FILE} .xyz)

XYZ_2_PDB ${INPUT_FILE}

CleanPDB ${INPUT_FILENAME}.pdb

AssignUniqueAtomsName ${INPUT_FILENAME}.pdb

ExtractMetal ${INPUT_FILENAME}_renamed.pdb

ExtractOrganicFragment ${INPUT_FILENAME}_renamed.pdb

PrepareMetalMol2 ${METAL}.pdb ${METAL_CHARGE} # comes from extractMetal()

PrepareOrganicFragment "${LIG_RESNAME}.pdb"

# Merge metal.pdb and lig.pdb keeping original atomID.
cat ${METAL}.pdb ${LIG_RESNAME}.pdb | sort -k2,2n > ${INPUT_FILENAME}_prep.pdb

LigatingAtoms ${INPUT_FILENAME}_prep.pdb ${METAL_ATOM_ID}
CheckFiles "PAIRS.txt"
PAIRS=$(<PAIRS.txt)

PrepareMCPBFile

echo "Done running prepare4mcpb.sh!"

#### End Main ####
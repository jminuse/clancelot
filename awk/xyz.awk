# This program removes all duplicate atoms in an xyz file under
# periodic boundary conditions, and outputs the correct number of atoms
# NOTE: The xyz file must be in fraction (crystal) coordinates!

# To use this script, type in command line:
# awk -f xyz.awk "old".xyz > "new".xyz

BEGIN {tol = 10e-4}    # criterion for deciding if two atoms overlaps
	
{
	if (NR>2 && NF==4) {
		for (rec=0; rec<=rec_new; rec++) {
			if (($2 - coord[rec+1 ":" 2])%1 >= -tol \
			&& ($2 - coord[rec+1 ":" 2])%1 <= tol \
			&& ($3 - coord[rec+1 ":" 3])%1 >= -tol \
			&& ($3 - coord[rec+1 ":" 3])%1 <= tol \
			&& ($4 - coord[rec+1 ":" 4])%1 >= -tol \
			&& ($4 - coord[rec+1 ":" 4])%1 <= tol)
				next;    # two atoms overlap
		}
		rec_new += 1;
		for (field=1; field<=NF; field++)
			coord[rec_new ":" field] = $field;
	}
}

END {
	printf("%d\n\n", rec_new);
	for (i=1; i<=rec_new; i++) {
		printf("%s\t%.4f\t%.4f\t%.4f\n", coord[i ":" 1], coord[i ":" 2],\
			coord[i ":" 3], coord[i ":" 4])
	}
}

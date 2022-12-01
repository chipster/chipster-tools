# TOOL split_tsv_file.R: "Split table" (Split TSV format table into several files.)
# INPUT table.tsv: "TSV file" TYPE GENERIC
# OUTPUT split{...}.tsv 
# PARAMETER substitute.header: "Substitute header" TYPE [yes, no] DEFAULT no (Substitute existing header with "identifier" and "sample")

if( substitute.header == "yes" ) {
	system("awk 'BEGIN{FS=\"\t\";}{if(NR==1){for(i=2;i<=NF;i++){gsub(/\r$/,\"\",$i);gsub(/\\s+$/,\"\",$i);gsub(/^\\s+/,\"\",$i);gsub(/\\.|\\/|\\\\/,\"\",$i);header[i]=$i;$i=\"sample\";};$1=\"identifier\"};for(i=2;i<=NF;i++){prefix=\"\";if(length(int((i-1)/n+1))==1){prefix=\"0\";};printf (i%n==0||i==7)?$1\"\t\"$i RS:$1\"\t\"$i FS > \"split\" prefix int((i-1)/n+1) \"_\" header[i] \".tsv\"}}' n=1 table.tsv")
} else {
	system("awk 'BEGIN{FS=\"\t\";}{for(i=2;i<=NF;i++){prefix=\"\";if(length(int((i-1)/n+1))==1){prefix=\"0\";};printf (i%n==0||i==7)?$1\"\t\"$i RS:$1\"\t\"$i FS > \"split\" prefix int((i-1)/n+1) \".tsv\"}}' n=1 table.tsv")
}


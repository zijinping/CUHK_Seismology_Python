# !/bin/sh
# Generate .list and .reverse from the header of the sac files on Old and New SCEDC format
# Created by jw, 2020/03/20

sac_folder='../tmp/EQ-sichuan/20190223213809'
#sac_folder='../tmp/EQ-sichuan/20190308111430'
output_name='sichuan'

# Example 1-3 Old SCEDC format and empty reverse file
saclst kstnm kcmpnm stla stlo stel knetwk f $sac_folder/*.SAC | \
awk '{ printf "%-4.4s %.3s %31.31s %9.4f %10.4f %5i %24.24s \n", $2, $3, " ", $4, $5, $6, $7 }'|sort -u > "${output_name}.list2"

# Example 4 New SCEDC format and empty reverse file
saclst knetwk kstnm kcmpnm stla stlo stel f $sac_folder/*.SAC | \
awk '{ printf "%2.2s  %5.5s %.33s %45.45s %9.4f %10.4f %5i \n", $2, $3, $4, " ", $5, $6, $7 }' |sort -u > "${output_name}.list4"
awk '{ printf "%-5.5s0        0\n", $2}' "${output_name}.list4" > "${output_name}.reverse"


# End message
echo Created: "${output_name}.list2" "${output_name}.list4"  "${output_name}.reverse"

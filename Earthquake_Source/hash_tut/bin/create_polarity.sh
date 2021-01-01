# !/bin/sh
# CREATE NEW SCEDC Phase format file 
# Created by jw, 2020/03/30


# Check arguments
if [ $# -ne 2 ]; then
echo Usage:\n $0 [station list] [outfile]
exit 0
fi

# User input station and polarity
echo station list: $1
echo Output file: $2 \n

echo "Station: [type q to quite]"
read sta

#Loop
while [ $sta != "q" ]; do
    echo Channels:
    grep $sta $1 | awk '{ print $3 }' 
    echo "Channel:"
    read cha
    echo Polarity:
    read polarity

    grep $sta $1 | grep $cha | awk '{ printf "%-5.5s%2.2s  %3.3s i %s \n", $2,$1,$3, polarity}' polarity="$polarity"
    grep $sta $1 | grep $cha | awk '{ printf "%-5.5s%2.2s  %3.3s i %s \n", $2,$1,$3, polarity}' polarity="$polarity" >>$2

    echo "\nStation: [type q to quite]"
    read sta
done

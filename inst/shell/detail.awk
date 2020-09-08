#/usr/bin/env awk
{
    print($1"\t"$2"\t"$3"\t"$4"\t"$7-$3"\t"$8-$3"\t"$9);
}
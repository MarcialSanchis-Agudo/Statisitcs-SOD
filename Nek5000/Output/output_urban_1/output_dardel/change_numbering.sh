a=1
for i in duct0.f*; do
  new=$(printf "duct0.f%05d" "$a") #04 pad to length of 4
  mv -i -- "$i" "$new"
  let a=a+1
done

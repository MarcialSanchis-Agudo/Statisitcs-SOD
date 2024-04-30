a=123
for i in animation*; do
  new=$(printf "animation_1.%04d.png" "$a") #04 pad to length of 4
  mv -i -- "$i" "$new"
  let a=a+1
done

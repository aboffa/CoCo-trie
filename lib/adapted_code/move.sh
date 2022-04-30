echo "moving bit_vector.hpp ..." &&
cp lib/adapted_code/bit_vector.hpp lib/ds2i/succinct &
echo "moving broadword.hpp ..." &&
cp lib/adapted_code/broadword.hpp lib/ds2i/succinct &
echo "moving intrinsics.hpp ..." &&
cp lib/adapted_code/intrinsics.hpp lib/ds2i/succinct &
echo "moving compact_elias_fano.hpp ..." &&
cp lib/adapted_code/compact_elias_fano.hpp lib/ds2i &
wait
echo "Done."

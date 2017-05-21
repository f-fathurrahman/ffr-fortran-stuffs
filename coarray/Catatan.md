# Catatan Coarray Fortran

Coarray Fortran adalah fitur khusus bahasa pemrograman Fortran
yang memungkinkan
programmer untuk menulis program paralel dalam Fortran.
Coarray Fortran didefinisikan dalam standard Fortran 2008.

Pada coarray Fortran, setiap proses disebut sebagai image.
Image ini dinomori dari 1, 2, 3, dan seterusnya.
(Hal ini berbeda dengan MPI, di mana proses dinomori dari 0, 1, 2, dst)

Setiap image mengeksekusi program yang sama. Karena kode program dapat bergantung
pada image mana yang mengeksekusinya, image-image dapat mengesekusi bagian yang
berbeda dari suatu program pada saat yang sama.



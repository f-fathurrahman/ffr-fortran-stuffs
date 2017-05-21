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

## Deklarasi coarray

Mendeklarasikan suatu variable, array atau objek sebagai coarray memungkinkan nilainya
pada suatu image dapat diakses pada image yang lain.

```fortran
REAL(8), CODIMENSION[*] :: ca(100)
! atau
REAL(8) :: ca(100)[*], cb(100)[*] ! coarray array-1d
```

ca dan cb adalah array-1d dengan 100 elemen. Kopi dari tiap array akan disimpan pada
setiap image. Attribut `CODIMENSION` menunjukkan bahwa nilai dari tiap array pada
tiap image dapat diakses oleh image lain.

```fortran
REAL(8) :: c2(0:9,4:12)[0:*]  ! coarray array-2d  
```

Pernyataan ini mendeklarasikan `c2` sebagai array-2d dan coarray-1d dengan batas bawah 0.
Batas atas pada kodimensi terakhir untuk tiap nonallocatable coarray harus `*`.

```fortran
INTEGER :: Npoints[*]  ! coarray skalar
```

Pada kode di atas, coarray `Npoints` adalah skalar dan nilainya pada tiap image dapat
diakses oleh image lain. Kodimensi dari coarray skalar selalu `*`.

Untuk array allocatable, semua dimensi dari kodimensi dari array harus `:`.

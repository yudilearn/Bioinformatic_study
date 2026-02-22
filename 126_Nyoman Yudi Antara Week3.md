Tugas 3

Analisis _Differential Expression Analysis_ (DEG)

Oleh: Nyoman Yudi Antara

**Pendahuluan**

Ekspresi gen merupakan proses biologis fundamental yang memungkinkan informasi genetik yang tersimpan dalam DNA diekspresikan menjadi molekul RNA dan selanjutnya diterjemahkan menjadi protein fungsional. Proses ini mencakup tahap transkripsi (DNA menjadi RNA) dan translasi (RNA menjadi protein), yang dikendalikan secara ketat oleh berbagai mekanisme regulasi seperti faktor transkripsi, modifikasi epigenetik, serta regulasi pascatranskripsi dan pascatranslasi. Regulasi ekspresi gen memastikan bahwa gen tertentu diaktifkan atau direpresi pada waktu, lokasi, dan tingkat yang sesuai dengan kebutuhan fisiologis sel. Meskipun seluruh sel dalam organisme multiseluler memiliki materi genetik yang sama, perbedaan pola ekspresi gen menghasilkan diferensiasi sel dan fungsi biologis yang beragam (NIH, 2026; Segundo-Val & Sanz-Lozano, 2016). Ekspresi gen dapat diketahui dengan beberapa metode seperti qPCR, Microarray, dan RNA sequencing yang ketiganya dapat menilai level atau kadar ekspresi suatu gen yang ingin diketahui. Studi komprehensif mengenai ekspresi gen disebut transkriptomik (Dong & Chen, 2013; Volgin, 2013).

Ekspresi gen atau _Differential Expression Analysis_ (DEG) dapat dilakukan menggunakan beberapa data seperti microarray ataupun RNA sequencing. Terdapat beberapa metode bioinformatik yang dapat digunakan untuk menganalisis ekspresi gen. Tugas ini ditujukan untuk mempelajari dan memahami proses analisis ekspresi gen atau _Differential Expression Analysis_ (DEG) menggunakan dataset yang ada di NCBI terutama GEO _(Gene Expression Omnibus)_.

**Metode**

Metode yang digunakan dalam tugas ini adalah menganalisis dataset yang ada di GEO mencari dan menggunakan dataset GEO dengan kode GSE 10072. GSE10072 merupakan dataset gene expression profiling dari pasien yang mengalami kanker paru-paru. Sampel terdiri dari dua kategori yaitu normal dan pasien yang mengalami kanker paru. Analisis tugas kali ini menggunakan bahasa pemrograman R sebagai platform untuk menganalisis differensial ekspresi gen yang dilanjutkan dengan menganalisis gene ontology (GO) dan analisis Pathway KEGG. Skrip hasil analisis selanjutnya dapat diakses di github: https://github.com/yudilearn/Bioinformatic_study.git

**Hasil**

Berdasarkan hasil analisis menghasilkan beberapa parameter sebagai berikut.

Gambar 1. Kurva distribusi nilai ekspresi

Gambar 1. Menunjukkan bahwa Kurva distribusi nilai ekspresi menunjukkan bahwa data telah ternormalisasi dengan sangat baik. Kedua kelompok (Adenocarcinoma dan Normal) memiliki puncak kepadatan (_density peak_) yang saling berhimpit pada rentang nilai Log2 yang sama. Hal ini mengonfirmasi bahwa perbedaan ekspresi yang ditemukan nantinya murni disebabkan oleh faktor biologis, bukan karena bias teknis atau perbedaan jumlah materi genetik antar sampel. Selanjutnya hasil analisis boxplot nilai ekspresi per sampel ditunjukkan pada Gambar 2.

B

A

Gambar 2. Distribusi nilai ekspresi per sampel (A), Reduksi Dimensi hasil ekspresi (B)

Berdasarkan Gambar 2A Analisis boxplot per sampel memperkuat hasil normalisasi tersebut. Median ekspresi dan rentang interkuartil (_interquartile range_) terlihat seragam di seluruh sampel (dari GSM254625 hingga GSM254729). Tidak ditemukan adanya sampel pencilan (_outlier_) yang signifikan, sehingga seluruh sampel layak untuk diikutsertakan dalam analisis diferensial ekspresi (DEG) lebih lanjut. Visualisasi menggunakan algoritma UMAP menunjukkan pemisahan klaster yang sangat kontras antara jaringan Adenocarcinoma (merah) dan jaringan Normal (biru). Sampel-sampel dalam kelompok yang sama berkumpul secara konsisten, yang mengindikasikan bahwa profil transkriptomik sel kanker paru sangat berbeda secara fundamental dibandingkan dengan jaringan paru yang sehat. Pemisahan yang tegas ini memberikan keyakinan tinggi terhadap validitas kelompok sampel yang diuji (Gambar 2B). Hasil ini dapat membantu melakukan analisis lebih lanjut yaitu melihat hasil ekspresi gen yang mengalami upregulated dan downregulated ditunjukkan dengan grafik vulcano plot dan heat map.

A B

Gambar 3. Gen yang mengalami upregulated dan downregulated (A) dan spesifik gen yang berbeda secara signifikan (B)

Gmbar 3A menunjukkan sejumlah besar gen yang mengalami perubahan ekspresi secara signifikan. Gen-gen yang berada di sisi kanan (merah) merupakan gen yang mengalami peningkatan ekspresi atau upregulated, sementara di sisi kiri (biru) adalah gen yang mengalami penurunan ekspresi atau downregulated pada penderita Adenocarcinoma. Dengan ambang batas log FC > 1 dan adj.P.Val < 0.01, terlihat sebaran gen yang masif, menunjukkan aktivitas molekuler yang sangat dinamis pada jaringan kanker. Terdapat 336 gen yang mengalami upregulated dan 584 gen yang mengalami downregulated (masing-masing 10 teratas ditampilkan pada tabel 1. Heatmap menampilkan pola ekspresi dari 50 gen teratas yang paling berbeda secara signifikan. Terlihat konsistensi pola di mana gen-gen seperti SPP1 menunjukkan penguatan ekspresi yang sangat tajam pada kelompok Adenocarcinoma (warna oranye/merah tua), sementara gen seperti AGER, CA4, dan FHL1 menunjukkan pelemahan ekspresi yang drastis dibandingkan jaringan normal (warna biru tua). Pengelompokan hierarki (_hierarchical clustering_) pada kolom juga berhasil mengelompokkan sampel Adenocarcinoma dan Normal secara sempurna tanpa kesalahan klasifikasi (Gambar 3B).

Tabel 1. Daftar 10 gen teratas yang mengalami upregulated dan down regulated

|     |     |     |     |     |
| --- | --- | --- | --- | --- |
| No  | Symbol | LogFc | adj.P.Val | Status gen |
| 1   | SPP1 | 4.364415 | 7.76E-37 | upregulated |
| 2   | COL11A1 | 3.061522 | 1.88E-21 | upregulated |
| 3   | MMP1 | 2.862036 | 5.80E-15 | upregulated |
| 4   | SPINK1 | 2.7584 | 1.35E-11 | upregulated |
| 5   | GREM1 | 2.548363 | 6.28E-23 | upregulated |
| 6   | COL10A1 | 2.463475 | 7.58E-20 | upregulated |
| 7   | TOP2A | 2.455963 | 2.43E-24 | upregulated |
| 8   | GREM1 | 2.441458 | 2.15E-21 | upregulated |
| 9   | TOX3 | 2.43796 | 3.82E-17 | upregulated |
| 10  | TOP2A | 2.422077 | 9.97E-22 | upregulated |
| 11  | AGER | \-4.41747 | 2.23E-37 | downregulated |
| 12  | SFTPC | \-4.04718 | 8.67E-24 | downregulated |
| 13  | SFTPC | \-4.00203 | 1.72E-23 | downregulated |
| 14  | SFTPC | \-3.99512 | 4.75E-24 | downregulated |
| 15  | SFTPC | \-3.98982 | 1.78E-33 | downregulated |
| 16  | FABP4 | \-3.83854 | 1.25E-35 | downregulated |
| 17  | CYP4B1 | \-3.7098 | 7.41E-24 | downregulated |
| 18  | WIF1 | \-3.68671 | 9.54E-24 | downregulated |
| 19  | ADH1B | \-3.67847 | 1.50E-26 | downregulated |
| 20  | TMEM100 | \-3.57031 | 1.79E-27 | downregulated |

Berdasarkan rangkaian analisis bioinformatika yang telah dilakukan, identifikasi gen yang terekspresi secara diferensial (_Differentially Expressed Genes_ atau DEG) antara kelompok Adenocarcinoma of the Lung dan Normal Lung Tissue menunjukkan perubahan transkriptomik yang masif. Hasil visualisasi melalui _Density Plot_ dan _Boxplot_ mengonfirmasi bahwa seluruh sampel telah ternormalisasi dengan baik, sehingga variasi ekspresi yang terdeteksi murni mencerminkan perbedaan kondisi patologis antar kelompok. Analisis reduksi dimensi menggunakan UMAP memberikan bukti kuat adanya diskriminasi profil genetik yang jelas, di mana sampel kanker dan normal terpisah ke dalam dua klaster yang sangat kontras. Hal ini divalidasi lebih lanjut melalui Volcano Plot yang memetakan distribusi gen berdasarkan signifikansi statistik ($adj.P.Val &lt; 0.01$) dan besarnya perubahan ekspresi ($\\log FC &gt; 1$). Dari pemetaan ini, ditemukan ratusan gen yang mengalami peningkatan ekspresi (_up-regulated_) dan penurunan ekspresi (_down-regulated_). Sebagai representasi visual dari gen yang paling berpengaruh, Heatmap Top 50 DEG menunjukkan pola ekspresi yang sangat konsisten di dalam masing-masing kelompok. Gen-gen seperti SPP1 muncul sebagai kandidat _up-regulated_ utama, sementara gen fungsional paru seperti AGER dan CA4 mengalami penurunan ekspresi yang signifikan pada jaringan kanker.

Namun, identifikasi daftar gen secara individu belum cukup untuk menjelaskan mekanisme patofisiologi kanker paru secara sistematis. Oleh karena itu, langkah analisis selanjutnya adalah melakukan Gene Ontology (GO) Enrichment Analysis untuk mengidentifikasi proses biologis yang paling dominan, serta KEGG Pathway Analysis untuk memetakan jalur signalisasi biokimia yang terganggu akibat perubahan ekspresi gen-gen tersebut (Rosati et al., 2024). Rangkaian analisis fungsional ini akan memberikan pemahaman yang lebih mendalam mengenai bagaimana regulasi gen tersebut berkontribusi terhadap progresi tumor, proliferasi sel, maupun kegagalan sistem pertahanan jaringan paru.

Gambar 4. Analisis GO dan KEGG

Hasil pengayaan Gene Ontology (GO) dan jalur KEGG pada gen yang ter-regulasi naik (_up-regulated_) mengungkapkan dominasi aktivitas mitosis dan gangguan siklus sel. Proses seperti _nuclear chromosome segregation_ dan jalur pensinyalan p53 menjadi indikator utama bahwa sel kanker telah kehilangan kendali atas stabilitas genomik dan terus melakukan proliferasi secara agresif. Fenomena ini menjelaskan pertumbuhan massa tumor yang cepat, di mana mesin pembelahan sel bekerja tanpa hambatan pos pemeriksaan (_checkpoint_) seluler yang normal (Chen et al., 2022). Sebaliknya, gen yang ter-regulasi turun (_down-regulated_) mencerminkan hancurnya integritas struktural dan pertahanan jaringan paru. Penurunan signifikan pada jalur _Focal adhesion, Adherens junction_, serta proses _regulation of angiogenesis_ menunjukkan bahwa sel kanker mulai kehilangan kontak dengan matriks ekstraseluler, yang menjadi prakondisi penting bagi invasi dan metastasis ke organ lain (Ben-Ze’ev, 1999). Selain itu, penekanan pada jalur _Complement and coagulation cascades_ mengindikasikan adanya strategi imunosupresi, di mana jaringan tumor menciptakan lingkungan yang memungkinkannya lolos dari pengawasan sistem kekebalan tubuh inang (Yang et al., 2024).

Referensi

Ben-Ze’ev, A. (1999). Focal Adhesions and Adherens Junctions: Their Role in Tumorigenesis. In E. E. Bittar, D. R. Garrod, A. J. North, & M. A. J. B. T.-A. in M. and C. B. Chidgey (Eds.), _The Adhesive Interaction of Cells_ (Vol. 28, pp. 135–163). Elsevier. https://doi.org/https://doi.org/10.1016/S1569-2558(08)60046-6

Chen, X., Zhang, T., Su, W., Dou, Z., Zhao, D., Jin, X., Lei, H., Wang, J., Xie, X., Cheng, B., Li, Q., Zhang, H., & Di, C. (2022). Mutant p53 in cancer: from molecular mechanism to therapeutic modulation. _Cell Death and Disease_, _13_(11). https://doi.org/10.1038/s41419-022-05408-1

Dong, Z. C., & Chen, Y. (2013). Transcriptomics: Advances and approaches. _Science China Life Sciences_, _56_(10), 960–967. https://doi.org/10.1007/s11427-013-4557-2

NIH. (2026). _Gene Expression_. National Human Genome Research Institute. https://www.genome.gov/genetics-glossary/Gene-Expression

Rosati, D., Palmieri, M., Brunelli, G., Morrione, A., Iannelli, F., Frullanti, E., & Giordano, A. (2024). Differential gene expression analysis pipelines and bioinformatic tools for the identification of specific biomarkers: A review. _Computational and Structural Biotechnology Journal_, _23_(February), 1154–1168. https://doi.org/10.1016/j.csbj.2024.02.018

Segundo-Val, I. S., & Sanz-Lozano, C. S. (2016). Introduction to the Gene Expression Analysis. _Methods in Molecular Biology (Clifton, N.J.)_, _1434_, 29–43. https://doi.org/10.1007/978-1-4939-3652-6_3

Volgin, D. V. (2013). Gene Expression: Analysis and Quantitation. _Animal Biotechnology: Models in Discovery and Translation_, 307–325. https://doi.org/10.1016/B978-0-12-416002-6.00017-1

Yang, J., Shen, L., Yang, J., Qu, Y., Gong, C., Zhou, F., Liu, Y., Luo, M., & Zhao, L. (2024). Complement and coagulation cascades are associated with prognosis and the immune microenvironment of lower-grade glioma. _Translational Cancer Research_, _13_(1), 112–136. https://doi.org/10.21037/tcr-23-906
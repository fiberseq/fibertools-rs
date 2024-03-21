# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and the project tries but probably doesn't to adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.4.2 (2024-03-21)

### Chore

 - <csr-id-d98cfa8cef7831d70fb4e8ebd2cccf38169e68cf/> clean docs
 - <csr-id-b4600cde832842ebdea80ca9485fb7f241d447e0/> clean docs
 - <csr-id-0222a23b3b412555d89af6e2dbc5a5160bc1bcf1/> clean changelog
 - <csr-id-fc4260d820084e0a9f321c2a681091b96770b93e/> update helps
 - <csr-id-fee819ae1b1d03d841abffd51a7698ae024443b1/> clean cargo config
 - <csr-id-069c6c95c575a60928998aab2e69044494ef3c55/> clean
 - <csr-id-a06a949d048a52600797d78fff8f35dfeccac136/> clean

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 9 commits contributed to the release over the course of 1 calendar day.
 - 1 day passed between releases.
 - 7 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 1 unique issue was worked on: [#47](https://github.com/fiberseq/fibertools-rs/issues/47)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#47](https://github.com/fiberseq/fibertools-rs/issues/47)**
    - 0.4.2 ([`bb7eff1`](https://github.com/fiberseq/fibertools-rs/commit/bb7eff193c7c47b290467b690cadde4066667487))
 * **Uncategorized**
    - Version ([`47be701`](https://github.com/fiberseq/fibertools-rs/commit/47be70198c732b506b8e66b42b77cde3080aa15b))
    - Clean docs ([`d98cfa8`](https://github.com/fiberseq/fibertools-rs/commit/d98cfa8cef7831d70fb4e8ebd2cccf38169e68cf))
    - Clean docs ([`b4600cd`](https://github.com/fiberseq/fibertools-rs/commit/b4600cde832842ebdea80ca9485fb7f241d447e0))
    - Clean changelog ([`0222a23`](https://github.com/fiberseq/fibertools-rs/commit/0222a23b3b412555d89af6e2dbc5a5160bc1bcf1))
    - Update helps ([`fc4260d`](https://github.com/fiberseq/fibertools-rs/commit/fc4260d820084e0a9f321c2a681091b96770b93e))
    - Clean cargo config ([`fee819a`](https://github.com/fiberseq/fibertools-rs/commit/fee819ae1b1d03d841abffd51a7698ae024443b1))
    - Clean ([`069c6c9`](https://github.com/fiberseq/fibertools-rs/commit/069c6c95c575a60928998aab2e69044494ef3c55))
    - Clean ([`a06a949`](https://github.com/fiberseq/fibertools-rs/commit/a06a949d048a52600797d78fff8f35dfeccac136))
</details>

## 0.4.1 (2024-03-20)

### Chore

### Chore

 - <csr-id-61967b45e2427eafdba6700c3aa4151ca4ba17c1/> make test data larger again since I remove the whole dir from cargo publish
 - <csr-id-93de92f6dbe3b8d5e596967f002d19ef0cd30458/> make test data larger again since I remove the whole dir from cargo publish
 - <csr-id-3a47fa8f33664c37c5ae9283b376360b7dff914d/> make test data larger again since I remove the whole dir from cargo publish
 - <csr-id-0a9ab3f816b0f52b8519b42b407577528dd14b19/> make test data larger again since I remove the whole dir from cargo publish
 - <csr-id-0b8f6a299120760452091ee3f14a4442d3876637/> make cargo publish smaller
 - <csr-id-2be96e24685037d23d7b9f1a8a4e18e53c5dda4a/> make cargo publish smaller
 - <csr-id-8a581721e0f45b6d6e774b38f7201142bd1f9eba/> make cargo publish smaller
 - <csr-id-1d7de118fbf56e6fbaa5e895f3ce0e1931c8f1be/> make cargo publish smaller
 - <csr-id-9a96cb31a28c25cbcd3f01a87b6f976071b0899e/> make cargo publish smaller
 - <csr-id-cbdda1a80fd3a06472b4ef3123949b5d1afa46ad/> make test data smaller
 - <csr-id-212b2c97b6e796aba879c92822640abe15561a36/> ignore notebook
 - <csr-id-285db7d2c044de140ed550bbe2bf617233b6f55e/> sync lock with branch
 - <csr-id-17b8cb6642bbee5f0e2813515d8256516dccb48d/> sync lock with branch

### New Features

 - <csr-id-547b8150f9159b41985dc59b77b254f148b200b0/> This release uses onnx files instead of pytorch files and moves the backend
   of m6a prediction to burn-rs instead of tch-rs.
   
   This allows us to use other ML backends in addition to pytorch, which
   means we can still run m6a prediction without the massive pytorch
   libraries and get the same results at the cost of some performance.
   
   I have also added tests for m6a predictions for each PacBio chemistry to
   ensure that the number and quality of predictions are identical to
   previous versions regardless of backend. Will also be a useful test in
   the future.
   
   This release also adds the new `ft-footprint` command to the `ft` CLI.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 24 commits contributed to the release over the course of 7 calendar days.
 - 20 days passed between releases.
 - 14 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 2 unique issues were worked on: [#43](https://github.com/fiberseq/fibertools-rs/issues/43), [#46](https://github.com/fiberseq/fibertools-rs/issues/46)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#43](https://github.com/fiberseq/fibertools-rs/issues/43)**
    - Footprinting ([`0429ff0`](https://github.com/fiberseq/fibertools-rs/commit/0429ff03a0f437f2d777d7b08b71efc930aa2898))
 * **[#46](https://github.com/fiberseq/fibertools-rs/issues/46)**
    - Convert pytorch models to onnx to allow more backends ([`547b815`](https://github.com/fiberseq/fibertools-rs/commit/547b8150f9159b41985dc59b77b254f148b200b0))
 * **Uncategorized**
    - Release fibertools-rs v0.4.1 ([`a123dac`](https://github.com/fiberseq/fibertools-rs/commit/a123dac23ba02d0d9d4ac8c4f8ffc4532f6d6516))
    - Make test data larger again since I remove the whole dir from cargo publish ([`61967b4`](https://github.com/fiberseq/fibertools-rs/commit/61967b45e2427eafdba6700c3aa4151ca4ba17c1))
    - Make test data larger again since I remove the whole dir from cargo publish ([`93de92f`](https://github.com/fiberseq/fibertools-rs/commit/93de92f6dbe3b8d5e596967f002d19ef0cd30458))
    - Make test data larger again since I remove the whole dir from cargo publish ([`3a47fa8`](https://github.com/fiberseq/fibertools-rs/commit/3a47fa8f33664c37c5ae9283b376360b7dff914d))
    - Make test data larger again since I remove the whole dir from cargo publish ([`0a9ab3f`](https://github.com/fiberseq/fibertools-rs/commit/0a9ab3f816b0f52b8519b42b407577528dd14b19))
    - Release fibertools-rs v0.4.0 ([`5f5fac5`](https://github.com/fiberseq/fibertools-rs/commit/5f5fac5c1336cba6d59dbeda464f34f8655e4070))
    - Release fibertools-rs v0.4.0 ([`01748de`](https://github.com/fiberseq/fibertools-rs/commit/01748de0fcaf07e67ecf330615389c3d29c7e9fe))
    - Make cargo publish smaller ([`0b8f6a2`](https://github.com/fiberseq/fibertools-rs/commit/0b8f6a299120760452091ee3f14a4442d3876637))
    - Release fibertools-rs v0.4.0 ([`3656e4b`](https://github.com/fiberseq/fibertools-rs/commit/3656e4b9ce3a0c3e7287f0abe04af72c608ecd98))
    - Release fibertools-rs v0.4.0 ([`e4815ad`](https://github.com/fiberseq/fibertools-rs/commit/e4815ad4c22acf94249063588ab48357fd2aab7b))
    - Make cargo publish smaller ([`2be96e2`](https://github.com/fiberseq/fibertools-rs/commit/2be96e24685037d23d7b9f1a8a4e18e53c5dda4a))
    - Release fibertools-rs v0.4.0 ([`4c505e0`](https://github.com/fiberseq/fibertools-rs/commit/4c505e0a165a527b2a7cb338d2414a263a11ce80))
    - Make cargo publish smaller ([`8a58172`](https://github.com/fiberseq/fibertools-rs/commit/8a581721e0f45b6d6e774b38f7201142bd1f9eba))
    - Release fibertools-rs v0.4.0 ([`31881bf`](https://github.com/fiberseq/fibertools-rs/commit/31881bfe94986035c343d762c3ecd67c051b2cbf))
    - Make cargo publish smaller ([`1d7de11`](https://github.com/fiberseq/fibertools-rs/commit/1d7de118fbf56e6fbaa5e895f3ce0e1931c8f1be))
    - Make cargo publish smaller ([`9a96cb3`](https://github.com/fiberseq/fibertools-rs/commit/9a96cb31a28c25cbcd3f01a87b6f976071b0899e))
    - Release fibertools-rs v0.4.0 ([`9e29925`](https://github.com/fiberseq/fibertools-rs/commit/9e2992599d89b2c71c8be31b711f471d8e08735a))
    - Make test data smaller ([`cbdda1a`](https://github.com/fiberseq/fibertools-rs/commit/cbdda1a80fd3a06472b4ef3123949b5d1afa46ad))
    - Release fibertools-rs v0.4.0 ([`5eba406`](https://github.com/fiberseq/fibertools-rs/commit/5eba4063dfa39103c712eb68ce583d906eb8c453))
    - Ignore notebook ([`212b2c9`](https://github.com/fiberseq/fibertools-rs/commit/212b2c97b6e796aba879c92822640abe15561a36))
    - Sync lock with branch ([`285db7d`](https://github.com/fiberseq/fibertools-rs/commit/285db7d2c044de140ed550bbe2bf617233b6f55e))
    - Sync lock with branch ([`17b8cb6`](https://github.com/fiberseq/fibertools-rs/commit/17b8cb6642bbee5f0e2813515d8256516dccb48d))
</details>

## 0.4.0 (2024-03-20)

0.4.0 failed to publish due to test data being too small. The release was re-published as 0.4.1 with the test data excluded from cargo publish.

## 0.3.9 (2024-02-28)

<csr-id-b0b9b6b8b11c7bd6be8906294897723f3b9c7167/>
<csr-id-4e524b1156a72d5787e626c5ddb9b884471d51a7/>
<csr-id-39133fdff9c3aca578c18385477fd8b21f7e4682/>
<csr-id-bbd76d25ab3f15cf5dcb7bebd1788dcc3a357815/>
<csr-id-f401efe4f3372fc6f1ee07daede818d73f9e9c77/>
<csr-id-24954a493bde458c43a2ab9cc1b1db92548ef43d/>

### Chore

 - <csr-id-b0b9b6b8b11c7bd6be8906294897723f3b9c7167/> Release fibertools-rs version 0.3.9
 - <csr-id-4e524b1156a72d5787e626c5ddb9b884471d51a7/> tch update
 - <csr-id-39133fdff9c3aca578c18385477fd8b21f7e4682/> tch update
 - <csr-id-bbd76d25ab3f15cf5dcb7bebd1788dcc3a357815/> pyft lock
 - <csr-id-f401efe4f3372fc6f1ee07daede818d73f9e9c77/> pyft lock
 - <csr-id-24954a493bde458c43a2ab9cc1b1db92548ef43d/> pyft lock

### New Features

 - <csr-id-6bd988f9daf0d751b8fbc3a1f844ddd5a1264e9b/> add optional min msp filter to fire output
 - <csr-id-c3a6e82a043938b3ffb83d1f2872030b516f6644/> add optional min msp filter to fire output
 - <csr-id-b0bbb5823899d95a8788d7632052368c737c05ed/> add optional min msp filter to fire output
 - <csr-id-32d8c57bdc92397b287ac430bce3783364c0a829/> add optional min msp filter to fire output

### Bug Fixes

 - <csr-id-77f5ad4405895499ec335f4d6193867bba3d1e98/> typo
 - <csr-id-f89c95be36c88dfeebdfe17d101573201ab192e0/> allow overwritting fire qual

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 12 commits contributed to the release over the course of 35 calendar days.
 - 35 days passed between releases.
 - 12 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.3.9 ([`b0b9b6b`](https://github.com/fiberseq/fibertools-rs/commit/b0b9b6b8b11c7bd6be8906294897723f3b9c7167))
    - Tch update ([`4e524b1`](https://github.com/fiberseq/fibertools-rs/commit/4e524b1156a72d5787e626c5ddb9b884471d51a7))
    - Tch update ([`39133fd`](https://github.com/fiberseq/fibertools-rs/commit/39133fdff9c3aca578c18385477fd8b21f7e4682))
    - Add optional min msp filter to fire output ([`6bd988f`](https://github.com/fiberseq/fibertools-rs/commit/6bd988f9daf0d751b8fbc3a1f844ddd5a1264e9b))
    - Add optional min msp filter to fire output ([`c3a6e82`](https://github.com/fiberseq/fibertools-rs/commit/c3a6e82a043938b3ffb83d1f2872030b516f6644))
    - Add optional min msp filter to fire output ([`b0bbb58`](https://github.com/fiberseq/fibertools-rs/commit/b0bbb5823899d95a8788d7632052368c737c05ed))
    - Add optional min msp filter to fire output ([`32d8c57`](https://github.com/fiberseq/fibertools-rs/commit/32d8c57bdc92397b287ac430bce3783364c0a829))
    - Typo ([`77f5ad4`](https://github.com/fiberseq/fibertools-rs/commit/77f5ad4405895499ec335f4d6193867bba3d1e98))
    - Allow overwritting fire qual ([`f89c95b`](https://github.com/fiberseq/fibertools-rs/commit/f89c95be36c88dfeebdfe17d101573201ab192e0))
    - Pyft lock ([`bbd76d2`](https://github.com/fiberseq/fibertools-rs/commit/bbd76d25ab3f15cf5dcb7bebd1788dcc3a357815))
    - Pyft lock ([`f401efe`](https://github.com/fiberseq/fibertools-rs/commit/f401efe4f3372fc6f1ee07daede818d73f9e9c77))
    - Pyft lock ([`24954a4`](https://github.com/fiberseq/fibertools-rs/commit/24954a493bde458c43a2ab9cc1b1db92548ef43d))
</details>

## 0.3.8 (2024-01-24)

<csr-id-81dcb9b0f0b43df614b7be981d72d04d8b9a3479/>
<csr-id-56b0db1305712fee0c17dc0d05077dc3f262a7d4/>
<csr-id-622a44591a713d644ec1cc9637e87e7cd7aa2262/>
<csr-id-b093f3aa4231ffc5c69b108e5032e86d2454ffaf/>
<csr-id-62f4f49fc3810d4d9d291b32ba067cc7af6de367/>
<csr-id-6393bd8692095aa77de54fb848d508769cda8480/>
<csr-id-dc7c74a837bb42cda668a22a4cd7e7589e88df08/>
<csr-id-616f8da1bc2b05d617faeaaf6d540e6d3b82e433/>
<csr-id-b6da9feab69869b0d75707b5067ac615d03dd64d/>

### Chore

 - <csr-id-81dcb9b0f0b43df614b7be981d72d04d8b9a3479/> pyft lock
 - <csr-id-56b0db1305712fee0c17dc0d05077dc3f262a7d4/> readme and help pages
 - <csr-id-622a44591a713d644ec1cc9637e87e7cd7aa2262/> readme and help pages
 - <csr-id-b093f3aa4231ffc5c69b108e5032e86d2454ffaf/> whitespace
 - <csr-id-62f4f49fc3810d4d9d291b32ba067cc7af6de367/> simplify colos
 - <csr-id-6393bd8692095aa77de54fb848d508769cda8480/> CI
 - <csr-id-dc7c74a837bb42cda668a22a4cd7e7589e88df08/> CI
 - <csr-id-616f8da1bc2b05d617faeaaf6d540e6d3b82e433/> clippy
 - <csr-id-b6da9feab69869b0d75707b5067ac615d03dd64d/> improved error msg

### New Features

 - <csr-id-9bac006224838961532c3754a9a24032860780c8/> center and extract now include fire scores in the range of 0-255
 - <csr-id-be89a6134cc2a56e7e08bcb3466736767bb33dac/> better fire paralization
 - <csr-id-ea06c18388ee4c60d36d88c1bb16ca8d03cc390d/> update pyft to have a writer function
 - <csr-id-c83b5a090dc7f24093bacea12ef638462da9251e/> update pyft to have a writer function
 - <csr-id-3bb2f45a437673cb70da5bfc3f115bdafc39f21c/> update pyft to have a writer function
 - <csr-id-e98ceed4e60f8bd75bda7c10eec174049b94b366/> update pyft

### Bug Fixes

 - <csr-id-16ef1b6c6047c5e3b963810280c1c62f2096d47c/> GBDT made a semvar breaking change, fixing and locking depandancies.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 20 commits contributed to the release over the course of 14 calendar days.
 - 14 days passed between releases.
 - 16 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 1 unique issue was worked on: [#38](https://github.com/fiberseq/fibertools-rs/issues/38)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#38](https://github.com/fiberseq/fibertools-rs/issues/38)**
    - Refactor ([`a24ce8e`](https://github.com/fiberseq/fibertools-rs/commit/a24ce8ee1b1b169371cfcf0fe9b378266045994d))
 * **Uncategorized**
    - Release fibertools-rs v0.3.8 ([`672bd9c`](https://github.com/fiberseq/fibertools-rs/commit/672bd9c70fcd43976ae53e7e751c96c6c3469c80))
    - Release fibertools-rs v0.3.7 ([`aaca6be`](https://github.com/fiberseq/fibertools-rs/commit/aaca6be84a656636ec87050a4cf6f38f6e7789da))
    - Pyft lock ([`81dcb9b`](https://github.com/fiberseq/fibertools-rs/commit/81dcb9b0f0b43df614b7be981d72d04d8b9a3479))
    - GBDT made a semvar breaking change, fixing and locking depandancies. ([`16ef1b6`](https://github.com/fiberseq/fibertools-rs/commit/16ef1b6c6047c5e3b963810280c1c62f2096d47c))
    - Readme and help pages ([`56b0db1`](https://github.com/fiberseq/fibertools-rs/commit/56b0db1305712fee0c17dc0d05077dc3f262a7d4))
    - Readme and help pages ([`622a445`](https://github.com/fiberseq/fibertools-rs/commit/622a44591a713d644ec1cc9637e87e7cd7aa2262))
    - Center and extract now include fire scores in the range of 0-255 ([`9bac006`](https://github.com/fiberseq/fibertools-rs/commit/9bac006224838961532c3754a9a24032860780c8))
    - Better fire paralization ([`be89a61`](https://github.com/fiberseq/fibertools-rs/commit/be89a6134cc2a56e7e08bcb3466736767bb33dac))
    - Whitespace ([`b093f3a`](https://github.com/fiberseq/fibertools-rs/commit/b093f3aa4231ffc5c69b108e5032e86d2454ffaf))
    - Simplify colos ([`62f4f49`](https://github.com/fiberseq/fibertools-rs/commit/62f4f49fc3810d4d9d291b32ba067cc7af6de367))
    - CI ([`6393bd8`](https://github.com/fiberseq/fibertools-rs/commit/6393bd8692095aa77de54fb848d508769cda8480))
    - CI ([`dc7c74a`](https://github.com/fiberseq/fibertools-rs/commit/dc7c74a837bb42cda668a22a4cd7e7589e88df08))
    - Clippy ([`616f8da`](https://github.com/fiberseq/fibertools-rs/commit/616f8da1bc2b05d617faeaaf6d540e6d3b82e433))
    - Update pyft to have a writer function ([`ea06c18`](https://github.com/fiberseq/fibertools-rs/commit/ea06c18388ee4c60d36d88c1bb16ca8d03cc390d))
    - Update pyft to have a writer function ([`c83b5a0`](https://github.com/fiberseq/fibertools-rs/commit/c83b5a090dc7f24093bacea12ef638462da9251e))
    - Update pyft to have a writer function ([`3bb2f45`](https://github.com/fiberseq/fibertools-rs/commit/3bb2f45a437673cb70da5bfc3f115bdafc39f21c))
    - Update pyft ([`e98ceed`](https://github.com/fiberseq/fibertools-rs/commit/e98ceed4e60f8bd75bda7c10eec174049b94b366))
    - Improved error msg ([`b6da9fe`](https://github.com/fiberseq/fibertools-rs/commit/b6da9feab69869b0d75707b5067ac615d03dd64d))
    - Unused ([`5824b74`](https://github.com/fiberseq/fibertools-rs/commit/5824b74d75673f95c39568f764bf05ee87ae9349))
</details>

## 0.3.7 (2024-01-10)

<csr-id-85b6a7180e8aa3f306deed66287ed8abd48c0cff/>
<csr-id-7b1e93362bce2d88c9af39fbe0e17bb2c71a01ba/>
<csr-id-76326aa133a1dfd84a14ef987f7275e3f9191476/>
<csr-id-1e08da2d9a5f5d4e6f8d26f46bc1d95c29f4692d/>
<csr-id-e6e45a3176369a21c8b7106c97a52644a2d596ce/>
<csr-id-cd1894501787e76c5590f471b43f7b6b441235ac/>
<csr-id-d33aed7543af326ca2096b5e675f9a83369741e5/>
<csr-id-df416375797b65a35055d8746083e522ece2fc57/>
<csr-id-ff5759efc4a11580d585abe52f4fbe97b44f948e/>
<csr-id-c7331897aaf9ec5b2a574aa07155675b7dee9a68/>
<csr-id-09249cbe1d46f4019a1943766b7a7c6ce22d732e/>
<csr-id-3cf7a30cd115133138bf097b56b81b728756b545/>
<csr-id-81dcb9b0f0b43df614b7be981d72d04d8b9a3479/>
<csr-id-56b0db1305712fee0c17dc0d05077dc3f262a7d4/>
<csr-id-622a44591a713d644ec1cc9637e87e7cd7aa2262/>
<csr-id-b093f3aa4231ffc5c69b108e5032e86d2454ffaf/>
<csr-id-62f4f49fc3810d4d9d291b32ba067cc7af6de367/>
<csr-id-6393bd8692095aa77de54fb848d508769cda8480/>
<csr-id-dc7c74a837bb42cda668a22a4cd7e7589e88df08/>
<csr-id-616f8da1bc2b05d617faeaaf6d540e6d3b82e433/>
<csr-id-b6da9feab69869b0d75707b5067ac615d03dd64d/>

### Chore

 - <csr-id-85b6a7180e8aa3f306deed66287ed8abd48c0cff/> add trace to ml mm parser
 - <csr-id-7b1e93362bce2d88c9af39fbe0e17bb2c71a01ba/> update.
 - <csr-id-76326aa133a1dfd84a14ef987f7275e3f9191476/> update.
 - <csr-id-1e08da2d9a5f5d4e6f8d26f46bc1d95c29f4692d/> update.
 - <csr-id-e6e45a3176369a21c8b7106c97a52644a2d596ce/> update.
 - <csr-id-cd1894501787e76c5590f471b43f7b6b441235ac/> update.
 - <csr-id-d33aed7543af326ca2096b5e675f9a83369741e5/> update.
 - <csr-id-df416375797b65a35055d8746083e522ece2fc57/> update.
 - <csr-id-ff5759efc4a11580d585abe52f4fbe97b44f948e/> update.

### Chore

 - <csr-id-81dcb9b0f0b43df614b7be981d72d04d8b9a3479/> pyft lock
 - <csr-id-56b0db1305712fee0c17dc0d05077dc3f262a7d4/> readme and help pages
 - <csr-id-622a44591a713d644ec1cc9637e87e7cd7aa2262/> readme and help pages
 - <csr-id-b093f3aa4231ffc5c69b108e5032e86d2454ffaf/> whitespace
 - <csr-id-62f4f49fc3810d4d9d291b32ba067cc7af6de367/> simplify colos
 - <csr-id-6393bd8692095aa77de54fb848d508769cda8480/> CI
 - <csr-id-dc7c74a837bb42cda668a22a4cd7e7589e88df08/> CI
 - <csr-id-616f8da1bc2b05d617faeaaf6d540e6d3b82e433/> clippy
 - <csr-id-b6da9feab69869b0d75707b5067ac615d03dd64d/> improved error msg

### New Features

 - <csr-id-7ddcf19be2814d616b3a337fcfe17944f0c5d1c5/> FIRE extract update
 - <csr-id-bfed25908b7380f535eaac05d53caa1402a3d097/> FIRE extract update
 - <csr-id-cebf54c21dac19a63e86d251b39b35e9e85ff7cb/> Adding a fire bed+ extract for the fire pipeline.
 - <csr-id-8db7b4187294429dbc1bb0eac02dcc9eeb14c9a7/> FIRE now works in fibertools! TODO add info to extract, etc.
 - <csr-id-3edee57927c8c92a883183da4ccfe26b7829f21e/> adding rle information to fire feats
 - <csr-id-b10a8c40689606c454554d4843322abafa292e91/> adding rle information to fire feats
 - <csr-id-c4591964da5251b8a93c41c6afa0bd368ffe08a6/> adding rle information to fire feats
 - <csr-id-3a96904f558f72245fdffe2817a0ca5a7b3afb3e/> adding rle information to fire feats
 - <csr-id-babb017266dcd3c0b6beb63d8f2da65f7b4c73fc/> adding rle information to fire feats
 - <csr-id-f830405d775e265790228400ac0c2fe5c2c33661/> adding rle information to fire feats
 - <csr-id-b954aab4022692cf5f9143bc8bb0cae067296f14/> adding rle information to fire feats
 - <csr-id-fd11b45c62f633c98872fcd3ca2250963eb75961/> adding rle information to fire feats
 - <csr-id-7c7c2c3a9f045f96f522ea407035b7d05d2000ae/> adding rle information to fire feats
 - <csr-id-e1f989722a1979d1fe1a0b3bde4fff1c2971c02d/> adding rle information to fire feats
 - <csr-id-7fac8c74bd322dce45b36abea71fffda64a70055/> adding rle information to fire feats
 - <csr-id-be21ec1fa6ab2de4fd89c9df9639c5c71a2d210a/> adding rle information to fire feats
 - <csr-id-20985ab8029090e9225f367faa5f920b32e085b7/> adding rle information to fire feats
 - <csr-id-0166d777a815aecd3d3d1e2f4e8d6e43c1172d8d/> adding rle information to fire feats
 - <csr-id-53ee8b40d44e89bec6f4e4c1fbc6013b326692f7/> adding rle information to fire feats
 - <csr-id-6ab4844c9e32c9a03337e36aa09f80462643c0cf/> adding rle information to fire feats
 - <csr-id-e4510e0e8ba22d74845d95531ccb185055447436/> adding rle information to fire feats
 - <csr-id-ee2732871e072271c9f18cb33458e59feb994da3/> adding rle information to fire feats
 - <csr-id-efb4372d2b59e82925a0c22d6002cd6e8e404a55/> adding rle information to fire feats
 - <csr-id-a7b98834cb8fdb0ecf6b832d11e90a0a1f492518/> update fire model
 - <csr-id-de02c2686f0ad3911123c392133ca680216a837d/> use threads better, fire feats
 - <csr-id-9ac73b76e0ee791b2f4445cff1f4c07b23ca2adf/> use threads better, fire feats
 - <csr-id-6d00af8092cb0a77fd59e5b7a7af5fd6356decbd/> use threads better, fire feats
 - <csr-id-fe239360d09012e96210ac16561a5493f72946e6/> use threads better, fire feats
 - <csr-id-16714115ebdf1c2ffaf7e0dd5d0a064996c8b186/> start working on decorators
 - <csr-id-1443446a9b676172735a9d6c29dc3d0e75ce7e81/> start working on decorators
 - <csr-id-9839e24972cd55c6786545824c5bd2b041e3ca7a/> FIRE io
 - <csr-id-51149cb14806a6a393b099ab845061f8ea3e7ee0/> Speed up writing FIRE results
 - <csr-id-bcca38ffdcc05c483387f81facbe2a7d10f064c7/> add a min ML options to add-nucs and keep working on FIRE predictions
 - <csr-id-1814ee48ed342a73920560a134b2a3bec5fd6c5e/> fire progress indicator.
 - <csr-id-77f60b934b3f1de1fae4d3de70af9db0d74e5a79/> Add ability to print FIRE features to a text file, TODO predict FIREs in rust.
 - <csr-id-564cc96d77ce10a40f8793afe6166dc2e397b5b3/> Move prediction into a feature only avalible through Cargo, since bioocnda wont let me do large installs anymore.
 - <csr-id-68550723237be3674783a86d11b3277ec09339d3/> add simplify options to center
 - <csr-id-c4d0c936f2d812d1f32c9455595f6141ff0b94fd/> add hp tag to center
 - <csr-id-eb5c54eb40d2b0efffa7fff30d607be4bfba9d66/> add hp tag to center
 - <csr-id-9bac006224838961532c3754a9a24032860780c8/> center and extract now include fire scores in the range of 0-255
 - <csr-id-be89a6134cc2a56e7e08bcb3466736767bb33dac/> better fire paralization
 - <csr-id-ea06c18388ee4c60d36d88c1bb16ca8d03cc390d/> update pyft to have a writer function
 - <csr-id-c83b5a090dc7f24093bacea12ef638462da9251e/> update pyft to have a writer function
 - <csr-id-3bb2f45a437673cb70da5bfc3f115bdafc39f21c/> update pyft to have a writer function
 - <csr-id-e98ceed4e60f8bd75bda7c10eec174049b94b366/> update pyft

### Bug Fixes

 - <csr-id-d50f3388aef6875b2f5a502874eb6d5712e61040/> allow for the unaligned read start and read end
 - <csr-id-51ea0b9661c9887e5e95375ef9af36ec04f1e371/> allow for the unaligned read start and read end
 - <csr-id-3a7134c519bee739b6258025d9e6888331f8494c/> make an iterator for fiberseq records
 - <csr-id-7ba5d14b07251c76071eb09eaa42cdf1df09e2d2/> try and better use threads.
 - <csr-id-43cc3dc7367b725f8e60f9d21c044dc638638f1f/> try and better use threads.
 - <csr-id-16ef1b6c6047c5e3b963810280c1c62f2096d47c/> GBDT made a semvar breaking change, fixing and locking depandancies.

### Other

 - <csr-id-c7331897aaf9ec5b2a574aa07155675b7dee9a68/> adding code outline for footprinting tool
 - <csr-id-09249cbe1d46f4019a1943766b7a7c6ce22d732e/> cli
 - <csr-id-3cf7a30cd115133138bf097b56b81b728756b545/> clippy

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 73 commits contributed to the release over the course of 81 calendar days.
 - 93 days passed between releases.
 - 56 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs v0.3.7 ([`9cdce73`](https://github.com/fiberseq/fibertools-rs/commit/9cdce73468536ec9d7720d517bd7cba70ea30d03))
    - Clippy ([`1f27489`](https://github.com/fiberseq/fibertools-rs/commit/1f2748914d155e8f55409df9ae9c161c5b3cf562))
    - Release fibertools-rs v0.3.6 ([`2f4a9ce`](https://github.com/fiberseq/fibertools-rs/commit/2f4a9ce4a5fc0c98fe78803eaad7aca3e5765d72))
    - Clippy ([`ca10f3e`](https://github.com/fiberseq/fibertools-rs/commit/ca10f3eb018e1186f4713e1ce47b2ba777c5ec63))
    - Paths ([`78e9b44`](https://github.com/fiberseq/fibertools-rs/commit/78e9b448b1ede8a14d02649e8e3360eed6521c3c))
    - Release fibertools-rs v0.3.6 ([`cf6891e`](https://github.com/fiberseq/fibertools-rs/commit/cf6891e61aaa2f14b382dc10e48214cb0c4e02ae))
    - FIRE extract update ([`7ddcf19`](https://github.com/fiberseq/fibertools-rs/commit/7ddcf19be2814d616b3a337fcfe17944f0c5d1c5))
    - FIRE extract update ([`bfed259`](https://github.com/fiberseq/fibertools-rs/commit/bfed25908b7380f535eaac05d53caa1402a3d097))
    - Adding code outline for footprinting tool ([`c733189`](https://github.com/fiberseq/fibertools-rs/commit/c7331897aaf9ec5b2a574aa07155675b7dee9a68))
    - Typo in wget ([`a7f24d8`](https://github.com/fiberseq/fibertools-rs/commit/a7f24d826f13a5d014e4de0952ac8fd2c4509cb1))
    - Update CI ([`6e3bf75`](https://github.com/fiberseq/fibertools-rs/commit/6e3bf754fedc5fc51b38405756b66417029aa307))
    - Add skip no m6a opts ([`a98c064`](https://github.com/fiberseq/fibertools-rs/commit/a98c06411dc19e443b7d6e519271150aebfff804))
    - Clippy ([`026862a`](https://github.com/fiberseq/fibertools-rs/commit/026862a80f3e2972cbe927adfca990ea741b2ad1))
    - Speed up extracting decorators, and move make decorators a struct ([`7176aa4`](https://github.com/fiberseq/fibertools-rs/commit/7176aa43e8631b2ea9ed4e63a912d7eca45a61af))
    - Adding a fire bed+ extract for the fire pipeline. ([`cebf54c`](https://github.com/fiberseq/fibertools-rs/commit/cebf54c21dac19a63e86d251b39b35e9e85ff7cb))
    - FIRE now works in fibertools! TODO add info to extract, etc. ([`8db7b41`](https://github.com/fiberseq/fibertools-rs/commit/8db7b4187294429dbc1bb0eac02dcc9eeb14c9a7))
    - Adding rle information to fire feats ([`3edee57`](https://github.com/fiberseq/fibertools-rs/commit/3edee57927c8c92a883183da4ccfe26b7829f21e))
    - Adding rle information to fire feats ([`b10a8c4`](https://github.com/fiberseq/fibertools-rs/commit/b10a8c40689606c454554d4843322abafa292e91))
    - Adding rle information to fire feats ([`c459196`](https://github.com/fiberseq/fibertools-rs/commit/c4591964da5251b8a93c41c6afa0bd368ffe08a6))
    - Adding rle information to fire feats ([`3a96904`](https://github.com/fiberseq/fibertools-rs/commit/3a96904f558f72245fdffe2817a0ca5a7b3afb3e))
    - Adding rle information to fire feats ([`babb017`](https://github.com/fiberseq/fibertools-rs/commit/babb017266dcd3c0b6beb63d8f2da65f7b4c73fc))
    - Adding rle information to fire feats ([`f830405`](https://github.com/fiberseq/fibertools-rs/commit/f830405d775e265790228400ac0c2fe5c2c33661))
    - Adding rle information to fire feats ([`b954aab`](https://github.com/fiberseq/fibertools-rs/commit/b954aab4022692cf5f9143bc8bb0cae067296f14))
    - Adding rle information to fire feats ([`fd11b45`](https://github.com/fiberseq/fibertools-rs/commit/fd11b45c62f633c98872fcd3ca2250963eb75961))
    - Adding rle information to fire feats ([`7c7c2c3`](https://github.com/fiberseq/fibertools-rs/commit/7c7c2c3a9f045f96f522ea407035b7d05d2000ae))
    - Adding rle information to fire feats ([`e1f9897`](https://github.com/fiberseq/fibertools-rs/commit/e1f989722a1979d1fe1a0b3bde4fff1c2971c02d))
    - Adding rle information to fire feats ([`7fac8c7`](https://github.com/fiberseq/fibertools-rs/commit/7fac8c74bd322dce45b36abea71fffda64a70055))
    - Adding rle information to fire feats ([`be21ec1`](https://github.com/fiberseq/fibertools-rs/commit/be21ec1fa6ab2de4fd89c9df9639c5c71a2d210a))
    - Adding rle information to fire feats ([`20985ab`](https://github.com/fiberseq/fibertools-rs/commit/20985ab8029090e9225f367faa5f920b32e085b7))
    - Adding rle information to fire feats ([`0166d77`](https://github.com/fiberseq/fibertools-rs/commit/0166d777a815aecd3d3d1e2f4e8d6e43c1172d8d))
    - Adding rle information to fire feats ([`53ee8b4`](https://github.com/fiberseq/fibertools-rs/commit/53ee8b40d44e89bec6f4e4c1fbc6013b326692f7))
    - Adding rle information to fire feats ([`6ab4844`](https://github.com/fiberseq/fibertools-rs/commit/6ab4844c9e32c9a03337e36aa09f80462643c0cf))
    - Adding rle information to fire feats ([`e4510e0`](https://github.com/fiberseq/fibertools-rs/commit/e4510e0e8ba22d74845d95531ccb185055447436))
    - Adding rle information to fire feats ([`ee27328`](https://github.com/fiberseq/fibertools-rs/commit/ee2732871e072271c9f18cb33458e59feb994da3))
    - Adding rle information to fire feats ([`efb4372`](https://github.com/fiberseq/fibertools-rs/commit/efb4372d2b59e82925a0c22d6002cd6e8e404a55))
    - Cli ([`09249cb`](https://github.com/fiberseq/fibertools-rs/commit/09249cbe1d46f4019a1943766b7a7c6ce22d732e))
    - Clippy ([`3cf7a30`](https://github.com/fiberseq/fibertools-rs/commit/3cf7a30cd115133138bf097b56b81b728756b545))
    - Update fire model ([`a7b9883`](https://github.com/fiberseq/fibertools-rs/commit/a7b98834cb8fdb0ecf6b832d11e90a0a1f492518))
    - Use threads better, fire feats ([`de02c26`](https://github.com/fiberseq/fibertools-rs/commit/de02c2686f0ad3911123c392133ca680216a837d))
    - Use threads better, fire feats ([`9ac73b7`](https://github.com/fiberseq/fibertools-rs/commit/9ac73b76e0ee791b2f4445cff1f4c07b23ca2adf))
    - Use threads better, fire feats ([`6d00af8`](https://github.com/fiberseq/fibertools-rs/commit/6d00af8092cb0a77fd59e5b7a7af5fd6356decbd))
    - Use threads better, fire feats ([`fe23936`](https://github.com/fiberseq/fibertools-rs/commit/fe239360d09012e96210ac16561a5493f72946e6))
    - Clean and move to FIRE ([`dea26a5`](https://github.com/fiberseq/fibertools-rs/commit/dea26a545075b26c00f52f1b62ee5ae815a62095))
    - Clean and move to FIRE ([`34b2c8f`](https://github.com/fiberseq/fibertools-rs/commit/34b2c8f795e21ae232850a03eb22c30a6aecea5d))
    - Clean and move to FIRE ([`1b3eccf`](https://github.com/fiberseq/fibertools-rs/commit/1b3eccfb30adf4c05c0d57981c7126f6a6273e73))
    - Clean and move to FIRE ([`77a90bb`](https://github.com/fiberseq/fibertools-rs/commit/77a90bbec650ca518c700250959939504f2c5fbf))
    - Start working on decorators ([`1671411`](https://github.com/fiberseq/fibertools-rs/commit/16714115ebdf1c2ffaf7e0dd5d0a064996c8b186))
    - Start working on decorators ([`1443446`](https://github.com/fiberseq/fibertools-rs/commit/1443446a9b676172735a9d6c29dc3d0e75ce7e81))
    - Add trace to ml mm parser ([`85b6a71`](https://github.com/fiberseq/fibertools-rs/commit/85b6a7180e8aa3f306deed66287ed8abd48c0cff))
    - FIRE io ([`9839e24`](https://github.com/fiberseq/fibertools-rs/commit/9839e24972cd55c6786545824c5bd2b041e3ca7a))
    - Speed up writing FIRE results ([`51149cb`](https://github.com/fiberseq/fibertools-rs/commit/51149cb14806a6a393b099ab845061f8ea3e7ee0))
    - Add a min ML options to add-nucs and keep working on FIRE predictions ([`bcca38f`](https://github.com/fiberseq/fibertools-rs/commit/bcca38ffdcc05c483387f81facbe2a7d10f064c7))
    - Allow for the unaligned read start and read end ([`d50f338`](https://github.com/fiberseq/fibertools-rs/commit/d50f3388aef6875b2f5a502874eb6d5712e61040))
    - Allow for the unaligned read start and read end ([`51ea0b9`](https://github.com/fiberseq/fibertools-rs/commit/51ea0b9661c9887e5e95375ef9af36ec04f1e371))
    - Make an iterator for fiberseq records ([`3a7134c`](https://github.com/fiberseq/fibertools-rs/commit/3a7134c519bee739b6258025d9e6888331f8494c))
    - Typo ([`f3af8c6`](https://github.com/fiberseq/fibertools-rs/commit/f3af8c671d929c8e9536eedcd35769f2af82db3c))
    - Try and better use threads. ([`7ba5d14`](https://github.com/fiberseq/fibertools-rs/commit/7ba5d14b07251c76071eb09eaa42cdf1df09e2d2))
    - Try and better use threads. ([`43cc3dc`](https://github.com/fiberseq/fibertools-rs/commit/43cc3dc7367b725f8e60f9d21c044dc638638f1f))
    - Fire progress indicator. ([`1814ee4`](https://github.com/fiberseq/fibertools-rs/commit/1814ee48ed342a73920560a134b2a3bec5fd6c5e))
    - Add ability to print FIRE features to a text file, TODO predict FIREs in rust. ([`77f60b9`](https://github.com/fiberseq/fibertools-rs/commit/77f60b934b3f1de1fae4d3de70af9db0d74e5a79))
    - Move prediction into a feature only avalible through Cargo, since bioocnda wont let me do large installs anymore. ([`564cc96`](https://github.com/fiberseq/fibertools-rs/commit/564cc96d77ce10a40f8793afe6166dc2e397b5b3))
    - Add simplify options to center ([`6855072`](https://github.com/fiberseq/fibertools-rs/commit/68550723237be3674783a86d11b3277ec09339d3))
    - Add hp tag to center ([`c4d0c93`](https://github.com/fiberseq/fibertools-rs/commit/c4d0c936f2d812d1f32c9455595f6141ff0b94fd))
    - Add hp tag to center ([`eb5c54e`](https://github.com/fiberseq/fibertools-rs/commit/eb5c54eb40d2b0efffa7fff30d607be4bfba9d66))
    - Start thinking about geting fire features ([`f163c02`](https://github.com/fiberseq/fibertools-rs/commit/f163c02277dd6267f6c2f7fba3574edb5a5cf86e))
    - Update. ([`7b1e933`](https://github.com/fiberseq/fibertools-rs/commit/7b1e93362bce2d88c9af39fbe0e17bb2c71a01ba))
    - Update. ([`76326aa`](https://github.com/fiberseq/fibertools-rs/commit/76326aa133a1dfd84a14ef987f7275e3f9191476))
    - Update. ([`1e08da2`](https://github.com/fiberseq/fibertools-rs/commit/1e08da2d9a5f5d4e6f8d26f46bc1d95c29f4692d))
    - Update. ([`e6e45a3`](https://github.com/fiberseq/fibertools-rs/commit/e6e45a3176369a21c8b7106c97a52644a2d596ce))
    - Update. ([`cd18945`](https://github.com/fiberseq/fibertools-rs/commit/cd1894501787e76c5590f471b43f7b6b441235ac))
    - Update. ([`d33aed7`](https://github.com/fiberseq/fibertools-rs/commit/d33aed7543af326ca2096b5e675f9a83369741e5))
    - Update. ([`df41637`](https://github.com/fiberseq/fibertools-rs/commit/df416375797b65a35055d8746083e522ece2fc57))
    - Update. ([`ff5759e`](https://github.com/fiberseq/fibertools-rs/commit/ff5759efc4a11580d585abe52f4fbe97b44f948e))
</details>

## 0.3.6 (2023-10-09)

<csr-id-32e9121f47e62e20b38600198f2f3fac3a4c88d6/>
<csr-id-c7331897aaf9ec5b2a574aa07155675b7dee9a68/>
<csr-id-09249cbe1d46f4019a1943766b7a7c6ce22d732e/>
<csr-id-3cf7a30cd115133138bf097b56b81b728756b545/>
<csr-id-85b6a7180e8aa3f306deed66287ed8abd48c0cff/>
<csr-id-7b1e93362bce2d88c9af39fbe0e17bb2c71a01ba/>
<csr-id-76326aa133a1dfd84a14ef987f7275e3f9191476/>
<csr-id-1e08da2d9a5f5d4e6f8d26f46bc1d95c29f4692d/>
<csr-id-e6e45a3176369a21c8b7106c97a52644a2d596ce/>
<csr-id-cd1894501787e76c5590f471b43f7b6b441235ac/>
<csr-id-d33aed7543af326ca2096b5e675f9a83369741e5/>
<csr-id-df416375797b65a35055d8746083e522ece2fc57/>
<csr-id-ff5759efc4a11580d585abe52f4fbe97b44f948e/>

### Chore

 - <csr-id-32e9121f47e62e20b38600198f2f3fac3a4c88d6/> update

### Other

 - <csr-id-c7331897aaf9ec5b2a574aa07155675b7dee9a68/> adding code outline for footprinting tool
 - <csr-id-09249cbe1d46f4019a1943766b7a7c6ce22d732e/> cli
 - <csr-id-3cf7a30cd115133138bf097b56b81b728756b545/> clippy

### Bug Fixes

 - <csr-id-d50f3388aef6875b2f5a502874eb6d5712e61040/> allow for the unaligned read start and read end
 - <csr-id-51ea0b9661c9887e5e95375ef9af36ec04f1e371/> allow for the unaligned read start and read end
 - <csr-id-3a7134c519bee739b6258025d9e6888331f8494c/> make an iterator for fiberseq records
 - <csr-id-7ba5d14b07251c76071eb09eaa42cdf1df09e2d2/> try and better use threads.
 - <csr-id-43cc3dc7367b725f8e60f9d21c044dc638638f1f/> try and better use threads.

### New Features

 - <csr-id-7ddcf19be2814d616b3a337fcfe17944f0c5d1c5/> FIRE extract update
 - <csr-id-bfed25908b7380f535eaac05d53caa1402a3d097/> FIRE extract update
 - <csr-id-cebf54c21dac19a63e86d251b39b35e9e85ff7cb/> Adding a fire bed+ extract for the fire pipeline.
 - <csr-id-8db7b4187294429dbc1bb0eac02dcc9eeb14c9a7/> FIRE now works in fibertools! TODO add info to extract, etc.
 - <csr-id-3edee57927c8c92a883183da4ccfe26b7829f21e/> adding rle information to fire feats
 - <csr-id-b10a8c40689606c454554d4843322abafa292e91/> adding rle information to fire feats
 - <csr-id-c4591964da5251b8a93c41c6afa0bd368ffe08a6/> adding rle information to fire feats
 - <csr-id-3a96904f558f72245fdffe2817a0ca5a7b3afb3e/> adding rle information to fire feats
 - <csr-id-babb017266dcd3c0b6beb63d8f2da65f7b4c73fc/> adding rle information to fire feats
 - <csr-id-f830405d775e265790228400ac0c2fe5c2c33661/> adding rle information to fire feats
 - <csr-id-b954aab4022692cf5f9143bc8bb0cae067296f14/> adding rle information to fire feats
 - <csr-id-fd11b45c62f633c98872fcd3ca2250963eb75961/> adding rle information to fire feats
 - <csr-id-7c7c2c3a9f045f96f522ea407035b7d05d2000ae/> adding rle information to fire feats
 - <csr-id-e1f989722a1979d1fe1a0b3bde4fff1c2971c02d/> adding rle information to fire feats
 - <csr-id-7fac8c74bd322dce45b36abea71fffda64a70055/> adding rle information to fire feats
 - <csr-id-be21ec1fa6ab2de4fd89c9df9639c5c71a2d210a/> adding rle information to fire feats
 - <csr-id-20985ab8029090e9225f367faa5f920b32e085b7/> adding rle information to fire feats
 - <csr-id-0166d777a815aecd3d3d1e2f4e8d6e43c1172d8d/> adding rle information to fire feats
 - <csr-id-53ee8b40d44e89bec6f4e4c1fbc6013b326692f7/> adding rle information to fire feats
 - <csr-id-6ab4844c9e32c9a03337e36aa09f80462643c0cf/> adding rle information to fire feats
 - <csr-id-e4510e0e8ba22d74845d95531ccb185055447436/> adding rle information to fire feats
 - <csr-id-ee2732871e072271c9f18cb33458e59feb994da3/> adding rle information to fire feats
 - <csr-id-efb4372d2b59e82925a0c22d6002cd6e8e404a55/> adding rle information to fire feats
 - <csr-id-a7b98834cb8fdb0ecf6b832d11e90a0a1f492518/> update fire model
 - <csr-id-de02c2686f0ad3911123c392133ca680216a837d/> use threads better, fire feats
 - <csr-id-9ac73b76e0ee791b2f4445cff1f4c07b23ca2adf/> use threads better, fire feats
 - <csr-id-6d00af8092cb0a77fd59e5b7a7af5fd6356decbd/> use threads better, fire feats
 - <csr-id-fe239360d09012e96210ac16561a5493f72946e6/> use threads better, fire feats
 - <csr-id-16714115ebdf1c2ffaf7e0dd5d0a064996c8b186/> start working on decorators
 - <csr-id-1443446a9b676172735a9d6c29dc3d0e75ce7e81/> start working on decorators
 - <csr-id-9839e24972cd55c6786545824c5bd2b041e3ca7a/> FIRE io
 - <csr-id-51149cb14806a6a393b099ab845061f8ea3e7ee0/> Speed up writing FIRE results
 - <csr-id-bcca38ffdcc05c483387f81facbe2a7d10f064c7/> add a min ML options to add-nucs and keep working on FIRE predictions
 - <csr-id-1814ee48ed342a73920560a134b2a3bec5fd6c5e/> fire progress indicator.
 - <csr-id-77f60b934b3f1de1fae4d3de70af9db0d74e5a79/> Add ability to print FIRE features to a text file, TODO predict FIREs in rust.
 - <csr-id-564cc96d77ce10a40f8793afe6166dc2e397b5b3/> Move prediction into a feature only avalible through Cargo, since bioocnda wont let me do large installs anymore.
 - <csr-id-68550723237be3674783a86d11b3277ec09339d3/> add simplify options to center
 - <csr-id-c4d0c936f2d812d1f32c9455595f6141ff0b94fd/> add hp tag to center
 - <csr-id-eb5c54eb40d2b0efffa7fff30d607be4bfba9d66/> add hp tag to center

### Chore

 - <csr-id-85b6a7180e8aa3f306deed66287ed8abd48c0cff/> add trace to ml mm parser
 - <csr-id-7b1e93362bce2d88c9af39fbe0e17bb2c71a01ba/> update.
 - <csr-id-76326aa133a1dfd84a14ef987f7275e3f9191476/> update.
 - <csr-id-1e08da2d9a5f5d4e6f8d26f46bc1d95c29f4692d/> update.
 - <csr-id-e6e45a3176369a21c8b7106c97a52644a2d596ce/> update.
 - <csr-id-cd1894501787e76c5590f471b43f7b6b441235ac/> update.
 - <csr-id-d33aed7543af326ca2096b5e675f9a83369741e5/> update.
 - <csr-id-df416375797b65a35055d8746083e522ece2fc57/> update.
 - <csr-id-ff5759efc4a11580d585abe52f4fbe97b44f948e/> update.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 2 commits contributed to the release.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs v0.3.6 ([`fd321d9`](https://github.com/fiberseq/fibertools-rs/commit/fd321d9a60f64f919061b68701e7c90a6f420140))
    - Update ([`32e9121`](https://github.com/fiberseq/fibertools-rs/commit/32e9121f47e62e20b38600198f2f3fac3a4c88d6))
</details>

## 0.3.5 (2023-10-09)

<csr-id-9b7ceba1f664fe6cdb553b5de43ff856a2cf94c2/>
<csr-id-c07f66b37f3e4e8d267578b7fc873c44e3df5d29/>

### Chore

 - <csr-id-9b7ceba1f664fe6cdb553b5de43ff856a2cf94c2/> improve progress bar
 - <csr-id-c07f66b37f3e4e8d267578b7fc873c44e3df5d29/> improve progress bar

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 6 commits contributed to the release.
 - 1 day passed between releases.
 - 2 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs v0.3.5 ([`1e351db`](https://github.com/fiberseq/fibertools-rs/commit/1e351dbacefc18429ab0b4ddc0f42d907cd1531d))
    - Release bio-io v0.3.2 ([`92b5918`](https://github.com/fiberseq/fibertools-rs/commit/92b5918cce1773a304017d4f5cee268145689580))
    - Improve progress bar ([`9b7ceba`](https://github.com/fiberseq/fibertools-rs/commit/9b7ceba1f664fe6cdb553b5de43ff856a2cf94c2))
    - Improve progress bar ([`c07f66b`](https://github.com/fiberseq/fibertools-rs/commit/c07f66b37f3e4e8d267578b7fc873c44e3df5d29))
    - Add notes ([`ff061a1`](https://github.com/fiberseq/fibertools-rs/commit/ff061a17a67ee8fa99e274613e4b2c4c1baa88a4))
    - Add notes ([`9fa61a0`](https://github.com/fiberseq/fibertools-rs/commit/9fa61a055952d16bd4b76010d5d3b44bf44d917b))
</details>

## 0.3.4 (2023-10-08)

<csr-id-e216e7d184ddcd235e0e22dd1da3e99382d78277/>

### Other

 - <csr-id-e216e7d184ddcd235e0e22dd1da3e99382d78277/> add docs.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 5 commits contributed to the release.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs v0.3.4 ([`8731ecb`](https://github.com/fiberseq/fibertools-rs/commit/8731ecb18abb4d1bcbcd4d95047c9f16785301f4))
    - Release bio-io v0.3.1 ([`efc9b3f`](https://github.com/fiberseq/fibertools-rs/commit/efc9b3f569490dac5b10462d359a834e324d62ef))
    - Release fibertools-rs v0.3.4, bio-io v0.3.1 ([`0697d8d`](https://github.com/fiberseq/fibertools-rs/commit/0697d8db1797a1a46bd276269be8687ebbd102df))
    - Add docs. ([`e216e7d`](https://github.com/fiberseq/fibertools-rs/commit/e216e7d184ddcd235e0e22dd1da3e99382d78277))
    - Release bio-io v0.3.1 ([`2d57bc1`](https://github.com/fiberseq/fibertools-rs/commit/2d57bc138cc952b466ac9b377169ba5994b31737))
</details>

## 0.3.3 (2023-10-08)

<csr-id-958e330daa4c51f7d69d9f7763333d78a92fb843/>
<csr-id-f60e8fa97ae89b8b19a52b5f9a0e62931262b622/>
<csr-id-d54dba8523a49eb9eac5d1b4743bcf9922094f8e/>
<csr-id-35d498ab1be9efa670d08c8191e639ebdbecbd5a/>
<csr-id-415cc7f2d151b18164e7cb2a8484e1369e61c2ca/>
<csr-id-e6ac4d3a7c63554a65511a15d4c7aca144203797/>
<csr-id-cef74228e5933e7b911db53b65352cfde5246b29/>
<csr-id-b1b148a865cbc0223d55ad2dc11b826f51210b54/>
<csr-id-8c818835c2c6f58103b3a8bfb0bbdd8e24379016/>
<csr-id-4ca7269d0851c70ea25c466bf2d8092075658e52/>
<csr-id-737fb50d4f311612f5f9aa5d1c342ed0d049b575/>
<csr-id-c6fec85f40be60b5b101ab7ed1ade439e07d4052/>
<csr-id-f2311eda7b3a400c39d93b8e95b7966b7ecba4da/>
<csr-id-5ffa686bc22216a497abfa48c5adf318a743a815/>
<csr-id-d166b29fdf08e24873afee8c37523293e288b762/>
<csr-id-30645d382bfc0a02f4ffa516044995a87acfea7c/>
<csr-id-897eca5bd467fc5ae7f888be1696af88bd3f4d9f/>
<csr-id-893f5f9e4dfca97b459281b5d667f7449a54364e/>

### Bug Fixes

 - <csr-id-6eb0527acb082d87eeb2dda2d898fe3642910dce/> progress bar display issues.

### Other

 - <csr-id-958e330daa4c51f7d69d9f7763333d78a92fb843/> add docs.
 - <csr-id-f60e8fa97ae89b8b19a52b5f9a0e62931262b622/> add docs.
 - <csr-id-d54dba8523a49eb9eac5d1b4743bcf9922094f8e/> add docs.
 - <csr-id-35d498ab1be9efa670d08c8191e639ebdbecbd5a/> add docs.
 - <csr-id-415cc7f2d151b18164e7cb2a8484e1369e61c2ca/> add docs.
 - <csr-id-e6ac4d3a7c63554a65511a15d4c7aca144203797/> add docs.
 - <csr-id-cef74228e5933e7b911db53b65352cfde5246b29/> add docs.
 - <csr-id-b1b148a865cbc0223d55ad2dc11b826f51210b54/> add docs.
 - <csr-id-8c818835c2c6f58103b3a8bfb0bbdd8e24379016/> add docs.
 - <csr-id-4ca7269d0851c70ea25c466bf2d8092075658e52/> add docs.
 - <csr-id-737fb50d4f311612f5f9aa5d1c342ed0d049b575/> add docs.
 - <csr-id-c6fec85f40be60b5b101ab7ed1ade439e07d4052/> add docs.
 - <csr-id-f2311eda7b3a400c39d93b8e95b7966b7ecba4da/> add docs.
 - <csr-id-5ffa686bc22216a497abfa48c5adf318a743a815/> add docs.
 - <csr-id-d166b29fdf08e24873afee8c37523293e288b762/> add docs.
 - <csr-id-30645d382bfc0a02f4ffa516044995a87acfea7c/> add docs.
 - <csr-id-897eca5bd467fc5ae7f888be1696af88bd3f4d9f/> add docs.
 - <csr-id-893f5f9e4dfca97b459281b5d667f7449a54364e/> add docs.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 20 commits contributed to the release.
 - 10 days passed between releases.
 - 19 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs v0.3.3 ([`85ee32a`](https://github.com/fiberseq/fibertools-rs/commit/85ee32aea21eccdfccc699ae1e2c4222e16855f7))
    - Add docs. ([`958e330`](https://github.com/fiberseq/fibertools-rs/commit/958e330daa4c51f7d69d9f7763333d78a92fb843))
    - Add docs. ([`f60e8fa`](https://github.com/fiberseq/fibertools-rs/commit/f60e8fa97ae89b8b19a52b5f9a0e62931262b622))
    - Add docs. ([`d54dba8`](https://github.com/fiberseq/fibertools-rs/commit/d54dba8523a49eb9eac5d1b4743bcf9922094f8e))
    - Add docs. ([`35d498a`](https://github.com/fiberseq/fibertools-rs/commit/35d498ab1be9efa670d08c8191e639ebdbecbd5a))
    - Add docs. ([`415cc7f`](https://github.com/fiberseq/fibertools-rs/commit/415cc7f2d151b18164e7cb2a8484e1369e61c2ca))
    - Add docs. ([`e6ac4d3`](https://github.com/fiberseq/fibertools-rs/commit/e6ac4d3a7c63554a65511a15d4c7aca144203797))
    - Add docs. ([`cef7422`](https://github.com/fiberseq/fibertools-rs/commit/cef74228e5933e7b911db53b65352cfde5246b29))
    - Add docs. ([`b1b148a`](https://github.com/fiberseq/fibertools-rs/commit/b1b148a865cbc0223d55ad2dc11b826f51210b54))
    - Add docs. ([`8c81883`](https://github.com/fiberseq/fibertools-rs/commit/8c818835c2c6f58103b3a8bfb0bbdd8e24379016))
    - Add docs. ([`4ca7269`](https://github.com/fiberseq/fibertools-rs/commit/4ca7269d0851c70ea25c466bf2d8092075658e52))
    - Add docs. ([`737fb50`](https://github.com/fiberseq/fibertools-rs/commit/737fb50d4f311612f5f9aa5d1c342ed0d049b575))
    - Add docs. ([`c6fec85`](https://github.com/fiberseq/fibertools-rs/commit/c6fec85f40be60b5b101ab7ed1ade439e07d4052))
    - Add docs. ([`f2311ed`](https://github.com/fiberseq/fibertools-rs/commit/f2311eda7b3a400c39d93b8e95b7966b7ecba4da))
    - Add docs. ([`5ffa686`](https://github.com/fiberseq/fibertools-rs/commit/5ffa686bc22216a497abfa48c5adf318a743a815))
    - Add docs. ([`d166b29`](https://github.com/fiberseq/fibertools-rs/commit/d166b29fdf08e24873afee8c37523293e288b762))
    - Add docs. ([`30645d3`](https://github.com/fiberseq/fibertools-rs/commit/30645d382bfc0a02f4ffa516044995a87acfea7c))
    - Add docs. ([`897eca5`](https://github.com/fiberseq/fibertools-rs/commit/897eca5bd467fc5ae7f888be1696af88bd3f4d9f))
    - Add docs. ([`893f5f9`](https://github.com/fiberseq/fibertools-rs/commit/893f5f9e4dfca97b459281b5d667f7449a54364e))
    - Progress bar display issues. ([`6eb0527`](https://github.com/fiberseq/fibertools-rs/commit/6eb0527acb082d87eeb2dda2d898fe3642910dce))
</details>

## 0.3.2 (2023-09-27)

<csr-id-12186d106224f211ac39844b6a8c82a54825ffcf/>
<csr-id-619081fa421d6658149f28b3a0458e83c4e16af2/>
<csr-id-fe3aa9a2511944b4078b3b7ed9dbb2ca756f9b81/>

### Chore

 - <csr-id-12186d106224f211ac39844b6a8c82a54825ffcf/> update deps
 - <csr-id-619081fa421d6658149f28b3a0458e83c4e16af2/> fix typo in cli message

### Bug Fixes

 - <csr-id-ea456be994b9209e32332bf5b46c183eba45ca3d/> extra comma in ft center wide format. fixes 21

### Other

 - <csr-id-fe3aa9a2511944b4078b3b7ed9dbb2ca756f9b81/> fix warning msg

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 8 commits contributed to the release over the course of 29 calendar days.
 - 45 days passed between releases.
 - 4 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs v0.3.2 ([`f39c11b`](https://github.com/fiberseq/fibertools-rs/commit/f39c11b54046af24684a8726ee7e6d54f52a6595))
    - Chore clippy ([`9c80dcb`](https://github.com/fiberseq/fibertools-rs/commit/9c80dcb2053ce42f2a33007954fc4506d50ce8a4))
    - Fix centered_query_start and centered_query_end values in ft center ([`a43c082`](https://github.com/fiberseq/fibertools-rs/commit/a43c0826c8a454e9ebb142e1fe4636a537883b6d))
    - Extra comma in ft center wide format. fixes 21 ([`ea456be`](https://github.com/fiberseq/fibertools-rs/commit/ea456be994b9209e32332bf5b46c183eba45ca3d))
    - Fix warning msg ([`fe3aa9a`](https://github.com/fiberseq/fibertools-rs/commit/fe3aa9a2511944b4078b3b7ed9dbb2ca756f9b81))
    - Update deps ([`12186d1`](https://github.com/fiberseq/fibertools-rs/commit/12186d106224f211ac39844b6a8c82a54825ffcf))
    - Fix typo in cli message ([`619081f`](https://github.com/fiberseq/fibertools-rs/commit/619081fa421d6658149f28b3a0458e83c4e16af2))
    - Add to cli message ([`7bf04d2`](https://github.com/fiberseq/fibertools-rs/commit/7bf04d280cef747391739896719b6fb0f84f6931))
</details>

## 0.3.1 (2023-08-13)

### Fixed

Added commas to wide format of ft center.

### New Features

 - <csr-id-c4c4b027dc400aeaf30b378270d073582e859a7e/> move progress bar into bamchunk iterator, unify progress bar.
 - <csr-id-b37eb2c6fcbe2b8e7c80e9ce380a24a951332887/> better progress bar for extract.
 - <csr-id-02070669179e6b6ff3383ade356d0a307e903905/> move progress bar into bamchunk iterator, unify progress bar.
 - <csr-id-d0a80ace111bb5e2477502ac264dd44fd988007e/> move progress bar into bamchunk iterator, unify progress bar.
 - <csr-id-678db90a111c54e33fed654031fe679b940d2d5a/> move progress bar into bamchunk iterator, unify progress bar.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 18 commits contributed to the release over the course of 13 calendar days.
 - 15 days passed between releases.
 - 5 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs v0.3.1 ([`6c5adf4`](https://github.com/fiberseq/fibertools-rs/commit/6c5adf4a88f04420e4024378cbed2bfdee73fe79))
    - Release fibertools-rs v0.3.1 ([`f5cc324`](https://github.com/fiberseq/fibertools-rs/commit/f5cc32460016b0a4ee9ea0804f236fc0aa73b803))
    - Update ([`9d870a9`](https://github.com/fiberseq/fibertools-rs/commit/9d870a928615c04cffe21b296970c9c6403d97de))
    - Release fibertools-rs v0.3.1 ([`0de730e`](https://github.com/fiberseq/fibertools-rs/commit/0de730e4949d61fff32fa06ff42c2d2f3bb16c29))
    - Update ([`1e7a6c1`](https://github.com/fiberseq/fibertools-rs/commit/1e7a6c1dfbb7805646a7d1e193b79e5699797269))
    - Release fibertools-rs v0.3.0 ([`b6cc95c`](https://github.com/fiberseq/fibertools-rs/commit/b6cc95c9b25b5da24fa522e2dd658f3a65d53bfe))
    - Release fibertools-rs v0.3.0 ([`3d3d534`](https://github.com/fiberseq/fibertools-rs/commit/3d3d534b7d214ae13a5b1e3f313007bbe21f385e))
    - Release fibertools-rs v0.3.0 ([`8262358`](https://github.com/fiberseq/fibertools-rs/commit/8262358293c29863d5eb44c6fe98df2c361659e2))
    - Update center.rs ([`90b1917`](https://github.com/fiberseq/fibertools-rs/commit/90b1917b38852a7bf3b3761cc08fa628a64e5a61))
    - Move progress bar into bamchunk iterator, unify progress bar. ([`0207066`](https://github.com/fiberseq/fibertools-rs/commit/02070669179e6b6ff3383ade356d0a307e903905))
    - Move progress bar into bamchunk iterator, unify progress bar. ([`d0a80ac`](https://github.com/fiberseq/fibertools-rs/commit/d0a80ace111bb5e2477502ac264dd44fd988007e))
    - Move progress bar into bamchunk iterator, unify progress bar. ([`678db90`](https://github.com/fiberseq/fibertools-rs/commit/678db90a111c54e33fed654031fe679b940d2d5a))
    - Move progress bar into bamchunk iterator, unify progress bar. ([`c4c4b02`](https://github.com/fiberseq/fibertools-rs/commit/c4c4b027dc400aeaf30b378270d073582e859a7e))
    - Better progress bar for extract. ([`b37eb2c`](https://github.com/fiberseq/fibertools-rs/commit/b37eb2c6fcbe2b8e7c80e9ce380a24a951332887))
    - Fmt ([`e50243a`](https://github.com/fiberseq/fibertools-rs/commit/e50243aa6f0857d188f9899203eda551bbdca6db))
    - Fmt ([`eccd328`](https://github.com/fiberseq/fibertools-rs/commit/eccd328d71a2b99bed146defe25ac2b0e1839d6c))
    - Fmt ([`250f080`](https://github.com/fiberseq/fibertools-rs/commit/250f080b65db4e4fd78f9edad6acf425a3ecdd78))
    - Fmt ([`9c10eff`](https://github.com/fiberseq/fibertools-rs/commit/9c10eff61ff5904ffcc85e214618b46b032ae616))
</details>

## 0.3.0 (2023-07-29)

<csr-id-6fe2beb91b2552292563fdae391c5ce026eded1d/>
<csr-id-df3ffba3aafa8403cdc95b6e0c31de07ac0a8114/>
<csr-id-96b4fc0638e991ec6c29dffbb1678b99ffc00b1d/>
<csr-id-1436dbfb05c542aa1f4c862d46113e43d0b396d8/>
<csr-id-04b1e9b02059718fea91b0aef71efc0d9d86f744/>
<csr-id-f73a36a716023835e16da50838752d74a2768843/>
<csr-id-168c9517b5137dbd58950f341d96a420c7d50935/>
<csr-id-9f711585d39c66ff45598d85dbf3783aca6ee44c/>
<csr-id-a80fe514ff3d4bf860dedbf7130a21d97b9537ab/>
<csr-id-c55fcc59749216e0f65c3d9c07efc14f32fb5fdf/>
<csr-id-2e1fe93f48b357fbb3459f5206fc5a994325d7a4/>
<csr-id-3e6aa945cb643e5d0559fac59b14a649e0de34b1/>
<csr-id-eea28d57a9eac79e3e3e0d9d16ede6dce51d242e/>
<csr-id-f5a65a612b8978f9742dab6990a8d52c6df7645a/>
<csr-id-b4278e4599a6634f37f05874a7e2e860c0b6a45a/>
<csr-id-6bb4c5c8b041c24743d2f89283c005f011084452/>
<csr-id-36a975a57afd5be38e4eb66429fc5b7f84ebe32d/>
<csr-id-15d535ae7dc90f86bec99a18015a22376c195788/>
<csr-id-aba1d7572a6109bb2e7d2acffea5f5e5c4f1ced6/>
<csr-id-9ea0b5f04f21cf74ab40a0ff1097eac996cc5b8b/>
<csr-id-69d87ea41febf0bd88380708cacea4232ae57f24/>
<csr-id-970c41f9ba4a8c26bf980f0ce87f635dd4f68bd3/>
<csr-id-49180bccd2dffdc49313e0b259dcdd0e80682bba/>
<csr-id-c8a2ab45cbff913c41bde48903c129b6e4d2f70a/>
<csr-id-8bf2168969acf67fa5069fb578f94941bd01e3f4/>

### Chore

- <csr-id-6fe2beb91b2552292563fdae391c5ce026eded1d/> simplifed center with new api.
- <csr-id-df3ffba3aafa8403cdc95b6e0c31de07ac0a8114/> refactor a large part of the code base to reduce redudance. TODO simplify center and basemods with this new api.
- <csr-id-96b4fc0638e991ec6c29dffbb1678b99ffc00b1d/> improve docs
- <csr-id-1436dbfb05c542aa1f4c862d46113e43d0b396d8/> drop darkmode
- <csr-id-04b1e9b02059718fea91b0aef71efc0d9d86f744/> fix docs and fetch
- <csr-id-f73a36a716023835e16da50838752d74a2768843/> regex update
- <csr-id-168c9517b5137dbd58950f341d96a420c7d50935/> readme
- <csr-id-9f711585d39c66ff45598d85dbf3783aca6ee44c/> clippy
- <csr-id-a80fe514ff3d4bf860dedbf7130a21d97b9537ab/> hide more deps under features for pyft
- <csr-id-c55fcc59749216e0f65c3d9c07efc14f32fb5fdf/> example py-ft
- <csr-id-2e1fe93f48b357fbb3459f5206fc5a994325d7a4/> unify api anming for liftover, optimize for inclusion in pyft
- <csr-id-3e6aa945cb643e5d0559fac59b14a649e0de34b1/> reduce number of copies for nuc and msp.
- <csr-id-eea28d57a9eac79e3e3e0d9d16ede6dce51d242e/> move nuc and msp logic out of extract into ranges struct and simplify. We ran on a whole genome ft extract to confirm results dont change.
- <csr-id-f5a65a612b8978f9742dab6990a8d52c6df7645a/> readme
- <csr-id-b4278e4599a6634f37f05874a7e2e860c0b6a45a/> clippy
- <csr-id-6bb4c5c8b041c24743d2f89283c005f011084452/> Release bio-io version 0.2.0
- <csr-id-36a975a57afd5be38e4eb66429fc5b7f84ebe32d/> Release bamlift version 0.2.0

### Chore

- <csr-id-8bf2168969acf67fa5069fb578f94941bd01e3f4/> readme

### Chore

- <csr-id-15d535ae7dc90f86bec99a18015a22376c195788/> readme
- <csr-id-aba1d7572a6109bb2e7d2acffea5f5e5c4f1ced6/> readme
- <csr-id-9ea0b5f04f21cf74ab40a0ff1097eac996cc5b8b/> readme
- <csr-id-69d87ea41febf0bd88380708cacea4232ae57f24/> readme
- <csr-id-970c41f9ba4a8c26bf980f0ce87f635dd4f68bd3/> readme
- <csr-id-49180bccd2dffdc49313e0b259dcdd0e80682bba/> readme
- <csr-id-c8a2ab45cbff913c41bde48903c129b6e4d2f70a/> readme

### New Features

 - <csr-id-b9f4950f2839c5d8df61593ec58997e439fc5c6a/> add qual to ft cetner, and clean the ft center code more
 - <csr-id-96f251b7f370f68e729357f52d3f2b6400d1af33/> adding ft center to the python module!
 - <csr-id-9da4fa3b27e991250185463de48ba0ebea866bce/> Reorganize pyft api and docs
 - <csr-id-248bf4e56c5b10c2dd4613888818a995e6f75529/> adding liftover to pyft
 - <csr-id-6a695ec7a71e5002709d5496c8aade541e737514/> very large improvment in speed of lifting over ranges, some liftover results differ from before by 1bp but it is rare.
 - <csr-id-9d07122e49b817f531b4a412ae77dced2e936f3a/> change the result of liftovers to be an Option<i64>.
 - <csr-id-02070669179e6b6ff3383ade356d0a307e903905/> move progress bar into bamchunk iterator, unify progress bar.
 - <csr-id-d0a80ace111bb5e2477502ac264dd44fd988007e/> move progress bar into bamchunk iterator, unify progress bar.
 - <csr-id-678db90a111c54e33fed654031fe679b940d2d5a/> move progress bar into bamchunk iterator, unify progress bar.
 - <csr-id-c4c4b027dc400aeaf30b378270d073582e859a7e/> move progress bar into bamchunk iterator, unify progress bar.
 - <csr-id-b37eb2c6fcbe2b8e7c80e9ce380a24a951332887/> better progress bar for extract.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 71 commits contributed to the release over the course of 7 calendar days.
 - 7 days passed between releases.
 - 31 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs v0.3.0 ([`740a4bc`](https://github.com/fiberseq/fibertools-rs/commit/740a4bce224e6345458450ea589c521c60fc998a))
    - Readme ([`d11301f`](https://github.com/fiberseq/fibertools-rs/commit/d11301fa6826563d397cb4ba7000810aae8f09c4))
    - Readme ([`b356b22`](https://github.com/fiberseq/fibertools-rs/commit/b356b22fa0596e2735c05cbf4004516e2899e411))
    - Readme ([`fe9dd43`](https://github.com/fiberseq/fibertools-rs/commit/fe9dd43268d0c5bc023146b14081346288d5d9be))
    - Readme ([`b532eb0`](https://github.com/fiberseq/fibertools-rs/commit/b532eb0d814cb965df7958f652279ecc8dd6f4fa))
    - Readme ([`2c4b6ef`](https://github.com/fiberseq/fibertools-rs/commit/2c4b6ef2eb7bb2d935e94be0f41ec727acf64c77))
    - Readme ([`ae12da9`](https://github.com/fiberseq/fibertools-rs/commit/ae12da9d5c7df82cfce89ed4e56d73200a7514e7))
    - Readme ([`9ca181d`](https://github.com/fiberseq/fibertools-rs/commit/9ca181de2ed1afb86e594b68a81118a50acd5b48))
    - Readme ([`5c2d004`](https://github.com/fiberseq/fibertools-rs/commit/5c2d004b326ae5a0164398d4bd436a1f2c6ad812))
    - Readme ([`f27c56d`](https://github.com/fiberseq/fibertools-rs/commit/f27c56d742a8eb8be2dd8d5271e5fa8568b56313))
    - Readme ([`67fc6eb`](https://github.com/fiberseq/fibertools-rs/commit/67fc6eb1d2eaa45d4e4fbd0fe750b1d4a22dd85e))
    - Release fibertools-rs v0.3.0 ([`668ff64`](https://github.com/fiberseq/fibertools-rs/commit/668ff648e125609e92d9c3b81d4f4b3481b8d152))
    - Release bio-io version 0.2.0 ([`6bb4c5c`](https://github.com/fiberseq/fibertools-rs/commit/6bb4c5c8b041c24743d2f89283c005f011084452))
    - Change ([`019e9de`](https://github.com/fiberseq/fibertools-rs/commit/019e9de0212ef414c69834b45c3dd3ebb29bcd5c))
    - Change ([`11c3915`](https://github.com/fiberseq/fibertools-rs/commit/11c39153fd5da5a5d630baca764b37eac884e872))
    - Change ([`20d9aa7`](https://github.com/fiberseq/fibertools-rs/commit/20d9aa70ff741fe198978abe4e8b8e46f7a9dd8e))
    - Change ([`de36de7`](https://github.com/fiberseq/fibertools-rs/commit/de36de75483774c8f89b5371c13d6000276a788b))
    - Release bamlift version 0.2.0 ([`36a975a`](https://github.com/fiberseq/fibertools-rs/commit/36a975a57afd5be38e4eb66429fc5b7f84ebe32d))
    - Change ([`70c9606`](https://github.com/fiberseq/fibertools-rs/commit/70c9606273dd05d1cafd7ef7e1a1ba49f1aa5aa5))
    - Change ([`49871ce`](https://github.com/fiberseq/fibertools-rs/commit/49871cef83c7229dd1754365876a7de33c1cdac1))
    - Change ([`cbf0de8`](https://github.com/fiberseq/fibertools-rs/commit/cbf0de8d74f0cb14b6d71b5bb6b3378666f29580))
    - Change ([`0a0c37b`](https://github.com/fiberseq/fibertools-rs/commit/0a0c37b4e62ed49048d78a836ddadcaa5911a844))
    - Change ([`4a5e5c4`](https://github.com/fiberseq/fibertools-rs/commit/4a5e5c4ecdcfcb82c86d9df0899ce3e3c4ef5e34))
    - Readme ([`99aeada`](https://github.com/fiberseq/fibertools-rs/commit/99aeadae8ba61da26d961a8f6318065e0a9ad442))
    - Improve fetch error msg ([`6332770`](https://github.com/fiberseq/fibertools-rs/commit/633277099e1fde8684733c05d6b192434ee01c10))
    - Update so pyft works with options results ([`a09faa2`](https://github.com/fiberseq/fibertools-rs/commit/a09faa2d68a6f13b2330b4a8e302fd40f3765c55))
    - Update so pyft works with options results ([`f225f9d`](https://github.com/fiberseq/fibertools-rs/commit/f225f9d8828b69b2af1a3c11b471fed84132e4d6))
    - Change the result of liftovers to be an Option<i64>. ([`9d07122`](https://github.com/fiberseq/fibertools-rs/commit/9d07122e49b817f531b4a412ae77dced2e936f3a))
    - Docs ([`07882ea`](https://github.com/fiberseq/fibertools-rs/commit/07882ea5396c178354ba735b651311a476c29cbc))
    - Docs ([`776713b`](https://github.com/fiberseq/fibertools-rs/commit/776713bf647ec6a361bb6bde34de412c6ec8af8c))
    - Change ([`ad0340e`](https://github.com/fiberseq/fibertools-rs/commit/ad0340e8d9f8b0ac744aea4785b9530c14572045))
    - Readme ([`6595c36`](https://github.com/fiberseq/fibertools-rs/commit/6595c3696dda1288f0f763fdc54d7311dfb790c9))
    - Add qual to ft cetner, and clean the ft center code more ([`b9f4950`](https://github.com/fiberseq/fibertools-rs/commit/b9f4950f2839c5d8df61593ec58997e439fc5c6a))
    - Docs ([`1a1f260`](https://github.com/fiberseq/fibertools-rs/commit/1a1f260323e262772ce7ecfaed2df52ab2597b5d))
    - Docs ([`2ab1645`](https://github.com/fiberseq/fibertools-rs/commit/2ab16455268ef311fbc60a17a7831e05b2f2cf8a))
    - Simplifed center with new api. ([`6fe2beb`](https://github.com/fiberseq/fibertools-rs/commit/6fe2beb91b2552292563fdae391c5ce026eded1d))
    - Refactor a large part of the code base to reduce redudance. TODO simplify center and basemods with this new api. ([`df3ffba`](https://github.com/fiberseq/fibertools-rs/commit/df3ffba3aafa8403cdc95b6e0c31de07ac0a8114))
    - Improve docs ([`96b4fc0`](https://github.com/fiberseq/fibertools-rs/commit/96b4fc0638e991ec6c29dffbb1678b99ffc00b1d))
    - Adding ft center to the python module! ([`96f251b`](https://github.com/fiberseq/fibertools-rs/commit/96f251b7f370f68e729357f52d3f2b6400d1af33))
    - Drop darkmode ([`1436dbf`](https://github.com/fiberseq/fibertools-rs/commit/1436dbfb05c542aa1f4c862d46113e43d0b396d8))
    - Fix docs and fetch ([`04b1e9b`](https://github.com/fiberseq/fibertools-rs/commit/04b1e9b02059718fea91b0aef71efc0d9d86f744))
    - Clean ([`1d5872d`](https://github.com/fiberseq/fibertools-rs/commit/1d5872dc510a370c065076e533af1913ccbc764f))
    - Reorganize pyft api and docs ([`9da4fa3`](https://github.com/fiberseq/fibertools-rs/commit/9da4fa3b27e991250185463de48ba0ebea866bce))
    - Check for sorted ([`8fbdc7c`](https://github.com/fiberseq/fibertools-rs/commit/8fbdc7c1a4af7039e62f249c49a4ca226023eefb))
    - Check for sorted ([`624dd33`](https://github.com/fiberseq/fibertools-rs/commit/624dd3317814594b47835be120a8d4bcadc775f6))
    - Block index ([`ab59d1d`](https://github.com/fiberseq/fibertools-rs/commit/ab59d1d3a650530c75d50b507a594ee7cc74eddc))
    - Regex update ([`f73a36a`](https://github.com/fiberseq/fibertools-rs/commit/f73a36a716023835e16da50838752d74a2768843))
    - Readme ([`8bf2168`](https://github.com/fiberseq/fibertools-rs/commit/8bf2168969acf67fa5069fb578f94941bd01e3f4))
    - Readme ([`15d535a`](https://github.com/fiberseq/fibertools-rs/commit/15d535ae7dc90f86bec99a18015a22376c195788))
    - Readme ([`aba1d75`](https://github.com/fiberseq/fibertools-rs/commit/aba1d7572a6109bb2e7d2acffea5f5e5c4f1ced6))
    - Readme ([`9ea0b5f`](https://github.com/fiberseq/fibertools-rs/commit/9ea0b5f04f21cf74ab40a0ff1097eac996cc5b8b))
    - Readme ([`69d87ea`](https://github.com/fiberseq/fibertools-rs/commit/69d87ea41febf0bd88380708cacea4232ae57f24))
    - Readme ([`970c41f`](https://github.com/fiberseq/fibertools-rs/commit/970c41f9ba4a8c26bf980f0ce87f635dd4f68bd3))
    - Readme ([`49180bc`](https://github.com/fiberseq/fibertools-rs/commit/49180bccd2dffdc49313e0b259dcdd0e80682bba))
    - Readme ([`c8a2ab4`](https://github.com/fiberseq/fibertools-rs/commit/c8a2ab45cbff913c41bde48903c129b6e4d2f70a))
    - Readme ([`168c951`](https://github.com/fiberseq/fibertools-rs/commit/168c9517b5137dbd58950f341d96a420c7d50935))
    - Readme ([`f5a65a6`](https://github.com/fiberseq/fibertools-rs/commit/f5a65a612b8978f9742dab6990a8d52c6df7645a))
    - Clippy ([`9f71158`](https://github.com/fiberseq/fibertools-rs/commit/9f711585d39c66ff45598d85dbf3783aca6ee44c))
    - Clippy ([`b4278e4`](https://github.com/fiberseq/fibertools-rs/commit/b4278e4599a6634f37f05874a7e2e860c0b6a45a))
    - Adding liftover to pyft ([`248bf4e`](https://github.com/fiberseq/fibertools-rs/commit/248bf4e56c5b10c2dd4613888818a995e6f75529))
    - Unify api anming for liftover, optimize for inclusion in pyft ([`2e1fe93`](https://github.com/fiberseq/fibertools-rs/commit/2e1fe93f48b357fbb3459f5206fc5a994325d7a4))
    - Unify api anming for liftover, optimize for inclusion in pyft ([`1c32392`](https://github.com/fiberseq/fibertools-rs/commit/1c32392c0bb24adfb3dc4bc78ab290e5405e2a62))
    - Clean ([`0e2a574`](https://github.com/fiberseq/fibertools-rs/commit/0e2a574cfd8e1b832a3f157b51c22a60713039a9))
    - Speed up and clean up approximate liftover and add a real test ([`c1ac44e`](https://github.com/fiberseq/fibertools-rs/commit/c1ac44e0a188a77900f9bf1dfef2dfbc42b3fdb9))
    - Update with ranges ([`cae1a41`](https://github.com/fiberseq/fibertools-rs/commit/cae1a415a2dac4aca1558d14a50eaa0cba258b09))
    - Very large improvment in speed of lifting over ranges, some liftover results differ from before by 1bp but it is rare. ([`6a695ec`](https://github.com/fiberseq/fibertools-rs/commit/6a695ec7a71e5002709d5496c8aade541e737514))
    - Hide more deps under features for pyft ([`a80fe51`](https://github.com/fiberseq/fibertools-rs/commit/a80fe514ff3d4bf860dedbf7130a21d97b9537ab))
    - Clean writes in center. ([`0d19ac1`](https://github.com/fiberseq/fibertools-rs/commit/0d19ac12281dd7fef0b571a745a3e6240bc771a3))
    - Reduce number of copies for nuc and msp. ([`3e6aa94`](https://github.com/fiberseq/fibertools-rs/commit/3e6aa945cb643e5d0559fac59b14a649e0de34b1))
    - Move nuc and msp logic out of extract into ranges struct and simplify. We ran on a whole genome ft extract to confirm results dont change. ([`eea28d5`](https://github.com/fiberseq/fibertools-rs/commit/eea28d57a9eac79e3e3e0d9d16ede6dce51d242e))
    - Example py-ft ([`c55fcc5`](https://github.com/fiberseq/fibertools-rs/commit/c55fcc59749216e0f65c3d9c07efc14f32fb5fdf))
</details>

## 0.2.6 (2023-07-21)

<csr-id-66b98f9f3c2644e79bf72de2ed070fba75b2ba38/>
<csr-id-9e8cda2fb47eff101eb3d458cf4599543e7d05a2/>
<csr-id-23e6b12a7601afb4aa21454048a53d803ca1bcdf/>
<csr-id-b03df77d99a8f28ecd6d43509a43ff92fb300329/>
<csr-id-6e8cd6f9f25dcabc15ef2a25595b3a3ee86de241/>
<csr-id-79ca845b0e72e832cdc57b1e1b58090a09417936/>
<csr-id-41bcd27e9de8654fd46bfb50d3c7efabad59a17e/>
<csr-id-66338fd576a66ee230aa4234d072a93b3a3b21d5/>
<csr-id-6f910c53a4e0cb7be411b724385e540517f21306/>
<csr-id-85d471df45a0c88567857427725852c491c041ff/>
<csr-id-c5eb60a8d42616eacfafb871cf4f46fb2fa924c8/>
<csr-id-5bd6803c346a409150daea7713878af46995e836/>
<csr-id-a59f41d0b2d4b97b277b738f04b23010be8d6449/>
<csr-id-7a633d71272a9425ae276d9b846fcda38f6526e5/>
<csr-id-8d16055edc4bf9bc5f94d5be42b86b2a2df712c8/>
<csr-id-467072ae6419b3ece8329f587467d0d0bf316911/>
<csr-id-09994bdbdd3949186d44ecef36199a8d228ea280/>
<csr-id-8ea9fb1efc4897e191468ee10aa66500456fa0fc/>
<csr-id-82399e6e702033acbb455ef6aea3a9935338d77f/>
<csr-id-fd115b7fdb6754d6f6885e901bc34d7364bfc7df/>
<csr-id-1d25caf7e8c19f97d39dec8eb69de7219e85621b/>
<csr-id-e73e00c849726e4577344ec8e964b3eda93ddbc6/>
<csr-id-a3f29cd798fcbeb071fce7ea7041c4e71a98f8f6/>
<csr-id-e85afd963f513fe578b2d0f3e368a7ae56b59d84/>
<csr-id-580968a9da3c2a66740d01dfcc0d8c574bfd46c8/>
<csr-id-6b1e714ccda049aef3199563c03d6f5ec1df9dfa/>
<csr-id-c8502ac7cf87d375a5aa9fcf3494a3cf8f43105c/>
<csr-id-b9274c2a7064cd615b079b3cc0b44d9e2f4b5e50/>
<csr-id-eb44bdaf8280689eb06e536266b64802f509a7b3/>
<csr-id-b98edf9f855ae24f81389a5ea481bf0b16dcd794/>
<csr-id-32f01826bbf06ac726464d59eee2769dfe80b138/>
<csr-id-94dfd521a3f637df47c996fd046d1de3cde85776/>
<csr-id-b3fa1c3b218379a3290fbed2382572b5e1c3e091/>
<csr-id-d4b75b1a569d694495bc1e87de651c1bb4e63778/>
<csr-id-73c918f9344fe362f602fa5356b11c9574c5ffec/>
<csr-id-5f635e04085c012d39d8afa362416ca3642fc7c5/>
<csr-id-c86cc585993506c55bd0243936767afe4a5f671c/>
<csr-id-15b825ed80a4b099c9ea9f21f64f1b695e20554a/>
<csr-id-26ef86389490b3161546271bea194e4c8a4a8b7f/>
<csr-id-d953400d9c7f676fea70fee51264c331ba7e1067/>
<csr-id-10bf0df63acca5812e3d2d3ad5134bc75abe73b4/>
<csr-id-cdde18e3b119ef36d57feed2f4740b7bde30c011/>
<csr-id-ac1987d33daf7fb7395ed1cf15b311390a15e60c/>
<csr-id-2f7c08e2b24a07054d547088ac2c96463a6229d3/>

### Bug Fixes

 - <csr-id-9fe3676efa853b79af85e2e1f6bb7849c013a0ec/> speed options
 - <csr-id-066ca9427211f5147f9528e66a342111e273296f/> rework fiberdata to consume bam record to avoid extra copy

### New Features

 - <csr-id-3ea32b01fda3e12f0cacae90d5a52495cec0b6bd/> add cpg to python modele and make easier access in the rust basemod api
 - <csr-id-1f6fea3cf07e1e798649df45d944a5396bf82902/> I have a minimal working python package yay!

### Chore

- <csr-id-b03df77d99a8f28ecd6d43509a43ff92fb300329/> rename iterator
- <csr-id-6e8cd6f9f25dcabc15ef2a25595b3a3ee86de241/> docs
- <csr-id-10bf0df63acca5812e3d2d3ad5134bc75abe73b4/> python docs
- <csr-id-cdde18e3b119ef36d57feed2f4740b7bde30c011/> readme
- <csr-id-ac1987d33daf7fb7395ed1cf15b311390a15e60c/> add to changelog

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 79 commits contributed to the release over the course of 2 calendar days.
 - 2 days passed between releases.
 - 48 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs v0.2.6 ([`3beecf6`](https://github.com/fiberseq/fibertools-rs/commit/3beecf697eafa9707abf5be30a2b773479125069))
    - Release fibertools-rs v0.2.6 ([`e973f17`](https://github.com/fiberseq/fibertools-rs/commit/e973f17def31e7c4731ae10143e6c5be5e01425a))
    - Release fibertools-rs v0.2.6 ([`ba8df9b`](https://github.com/fiberseq/fibertools-rs/commit/ba8df9b8e918be85beef13411cb3aab588634056))
    - Release fibertools-rs v0.2.6 ([`1b5dcfa`](https://github.com/fiberseq/fibertools-rs/commit/1b5dcfa18eef69a8e4a3f7e3e514bb5d5ae87b7f))
    - Release fibertools-rs v0.2.6 ([`5aa1908`](https://github.com/fiberseq/fibertools-rs/commit/5aa19089b0cf0614a7a8d0b4f8dc477d6e808ea5))
    - Release fibertools-rs v0.2.6 ([`656682a`](https://github.com/fiberseq/fibertools-rs/commit/656682a1b70871df4ca52484351d8296755cbe02))
    - Release fibertools-rs v0.2.6 ([`2d2f0f9`](https://github.com/fiberseq/fibertools-rs/commit/2d2f0f96800779088138f6671ab81acfb61cdacc))
    - Release fibertools-rs v0.2.6 ([`b80d856`](https://github.com/fiberseq/fibertools-rs/commit/b80d8567ee7ee899734fa5d2c5fd2b6604d45c78))
    - Release fibertools-rs v0.2.6 ([`14562fd`](https://github.com/fiberseq/fibertools-rs/commit/14562fd4b9552752ddba160cd3cba04e5f90102b))
    - Try again ([`4aa06eb`](https://github.com/fiberseq/fibertools-rs/commit/4aa06eb00c66afb2e5c0273247b81dbabfa15417))
    - Release fibertools-rs v0.2.6 ([`e247747`](https://github.com/fiberseq/fibertools-rs/commit/e24774703e6addefabcbe4b1a4df0a91383bdfa8))
    - Try again ([`09da733`](https://github.com/fiberseq/fibertools-rs/commit/09da733566d266063c7d74cbd610da2fe33b5f2f))
    - Release fibertools-rs v0.2.6 ([`9b1c3f2`](https://github.com/fiberseq/fibertools-rs/commit/9b1c3f2241dce70ece8ded8feeed19251687d1b1))
    - Release fibertools-rs v0.2.6 ([`0c525b5`](https://github.com/fiberseq/fibertools-rs/commit/0c525b5e0d6d81a1a82dd14a7fb6ea3c500c8c35))
    - Py-ft release ([`4de2b11`](https://github.com/fiberseq/fibertools-rs/commit/4de2b11dbcdec8836b53d0abef5444f87fd8c372))
    - Py-ft release ([`b909bc8`](https://github.com/fiberseq/fibertools-rs/commit/b909bc8df6adda1e0cef2ce6198c14e44ba6cbeb))
    - Larger chunks not worth ([`911bf40`](https://github.com/fiberseq/fibertools-rs/commit/911bf40a22eaace2afe659efc1000808078792fa))
    - Larger chunks ([`171cf51`](https://github.com/fiberseq/fibertools-rs/commit/171cf5104464aea984c8c1266835a2f385e49fa0))
    - Speed options ([`9fe3676`](https://github.com/fiberseq/fibertools-rs/commit/9fe3676efa853b79af85e2e1f6bb7849c013a0ec))
    - Rework fiberdata to consume bam record to avoid extra copy ([`066ca94`](https://github.com/fiberseq/fibertools-rs/commit/066ca9427211f5147f9528e66a342111e273296f))
    - Clippy ([`ecd0a61`](https://github.com/fiberseq/fibertools-rs/commit/ecd0a61bf1518476fa19be4c77554c68076b96fd))
    - Update rusthtslib version ([`d262afa`](https://github.com/fiberseq/fibertools-rs/commit/d262afac986cf53e97f3a53d41c87d605231c03e))
    - Move bam tags into bio_io ([`abcbcda`](https://github.com/fiberseq/fibertools-rs/commit/abcbcda65a9006c21c892bfdff509cf6e7e55c55))
    - Add logging and iterator and aligned blocks ([`7671998`](https://github.com/fiberseq/fibertools-rs/commit/76719983d8104a04cabd9846a5529f8e65fc7f96))
    - Rename iterator ([`66b98f9`](https://github.com/fiberseq/fibertools-rs/commit/66b98f9f3c2644e79bf72de2ed070fba75b2ba38))
    - Rename iterator ([`9e8cda2`](https://github.com/fiberseq/fibertools-rs/commit/9e8cda2fb47eff101eb3d458cf4599543e7d05a2))
    - Rename iterator ([`23e6b12`](https://github.com/fiberseq/fibertools-rs/commit/23e6b12a7601afb4aa21454048a53d803ca1bcdf))
    - Rename iterator ([`b03df77`](https://github.com/fiberseq/fibertools-rs/commit/b03df77d99a8f28ecd6d43509a43ff92fb300329))
    - Docs ([`6e8cd6f`](https://github.com/fiberseq/fibertools-rs/commit/6e8cd6f9f25dcabc15ef2a25595b3a3ee86de241))
    - Docs ([`79ca845`](https://github.com/fiberseq/fibertools-rs/commit/79ca845b0e72e832cdc57b1e1b58090a09417936))
    - Add cpg to python modele and make easier access in the rust basemod api ([`3ea32b0`](https://github.com/fiberseq/fibertools-rs/commit/3ea32b01fda3e12f0cacae90d5a52495cec0b6bd))
    - Docs ([`41bcd27`](https://github.com/fiberseq/fibertools-rs/commit/41bcd27e9de8654fd46bfb50d3c7efabad59a17e))
    - Docs ([`66338fd`](https://github.com/fiberseq/fibertools-rs/commit/66338fd576a66ee230aa4234d072a93b3a3b21d5))
    - Docs ([`6f910c5`](https://github.com/fiberseq/fibertools-rs/commit/6f910c53a4e0cb7be411b724385e540517f21306))
    - Docs ([`85d471d`](https://github.com/fiberseq/fibertools-rs/commit/85d471df45a0c88567857427725852c491c041ff))
    - Docs ([`c5eb60a`](https://github.com/fiberseq/fibertools-rs/commit/c5eb60a8d42616eacfafb871cf4f46fb2fa924c8))
    - Docs ([`5bd6803`](https://github.com/fiberseq/fibertools-rs/commit/5bd6803c346a409150daea7713878af46995e836))
    - Docs ([`a59f41d`](https://github.com/fiberseq/fibertools-rs/commit/a59f41d0b2d4b97b277b738f04b23010be8d6449))
    - Docs ([`7a633d7`](https://github.com/fiberseq/fibertools-rs/commit/7a633d71272a9425ae276d9b846fcda38f6526e5))
    - Docs ([`8d16055`](https://github.com/fiberseq/fibertools-rs/commit/8d16055edc4bf9bc5f94d5be42b86b2a2df712c8))
    - Docs ([`467072a`](https://github.com/fiberseq/fibertools-rs/commit/467072ae6419b3ece8329f587467d0d0bf316911))
    - Docs ([`09994bd`](https://github.com/fiberseq/fibertools-rs/commit/09994bdbdd3949186d44ecef36199a8d228ea280))
    - Docs ([`8ea9fb1`](https://github.com/fiberseq/fibertools-rs/commit/8ea9fb1efc4897e191468ee10aa66500456fa0fc))
    - Docs ([`82399e6`](https://github.com/fiberseq/fibertools-rs/commit/82399e6e702033acbb455ef6aea3a9935338d77f))
    - Docs ([`fd115b7`](https://github.com/fiberseq/fibertools-rs/commit/fd115b7fdb6754d6f6885e901bc34d7364bfc7df))
    - Docs ([`1d25caf`](https://github.com/fiberseq/fibertools-rs/commit/1d25caf7e8c19f97d39dec8eb69de7219e85621b))
    - Docs ([`e73e00c`](https://github.com/fiberseq/fibertools-rs/commit/e73e00c849726e4577344ec8e964b3eda93ddbc6))
    - Docs ([`a3f29cd`](https://github.com/fiberseq/fibertools-rs/commit/a3f29cd798fcbeb071fce7ea7041c4e71a98f8f6))
    - Docs ([`e85afd9`](https://github.com/fiberseq/fibertools-rs/commit/e85afd963f513fe578b2d0f3e368a7ae56b59d84))
    - Docs ([`580968a`](https://github.com/fiberseq/fibertools-rs/commit/580968a9da3c2a66740d01dfcc0d8c574bfd46c8))
    - Docs ([`6b1e714`](https://github.com/fiberseq/fibertools-rs/commit/6b1e714ccda049aef3199563c03d6f5ec1df9dfa))
    - Docs ([`c8502ac`](https://github.com/fiberseq/fibertools-rs/commit/c8502ac7cf87d375a5aa9fcf3494a3cf8f43105c))
    - Docs ([`b9274c2`](https://github.com/fiberseq/fibertools-rs/commit/b9274c2a7064cd615b079b3cc0b44d9e2f4b5e50))
    - Docs ([`eb44bda`](https://github.com/fiberseq/fibertools-rs/commit/eb44bdaf8280689eb06e536266b64802f509a7b3))
    - Rtd theme ([`44f93a4`](https://github.com/fiberseq/fibertools-rs/commit/44f93a42b1d5cb30fefad99e1b4dcf19f4f7f0c0))
    - Python docs ([`b98edf9`](https://github.com/fiberseq/fibertools-rs/commit/b98edf9f855ae24f81389a5ea481bf0b16dcd794))
    - Python docs ([`32f0182`](https://github.com/fiberseq/fibertools-rs/commit/32f01826bbf06ac726464d59eee2769dfe80b138))
    - Python docs ([`94dfd52`](https://github.com/fiberseq/fibertools-rs/commit/94dfd521a3f637df47c996fd046d1de3cde85776))
    - Python docs ([`b3fa1c3`](https://github.com/fiberseq/fibertools-rs/commit/b3fa1c3b218379a3290fbed2382572b5e1c3e091))
    - Python docs ([`d4b75b1`](https://github.com/fiberseq/fibertools-rs/commit/d4b75b1a569d694495bc1e87de651c1bb4e63778))
    - Python docs ([`73c918f`](https://github.com/fiberseq/fibertools-rs/commit/73c918f9344fe362f602fa5356b11c9574c5ffec))
    - Python docs ([`5f635e0`](https://github.com/fiberseq/fibertools-rs/commit/5f635e04085c012d39d8afa362416ca3642fc7c5))
    - Python docs ([`c86cc58`](https://github.com/fiberseq/fibertools-rs/commit/c86cc585993506c55bd0243936767afe4a5f671c))
    - Python docs ([`15b825e`](https://github.com/fiberseq/fibertools-rs/commit/15b825ed80a4b099c9ea9f21f64f1b695e20554a))
    - Python docs ([`26ef863`](https://github.com/fiberseq/fibertools-rs/commit/26ef86389490b3161546271bea194e4c8a4a8b7f))
    - Python docs ([`d953400`](https://github.com/fiberseq/fibertools-rs/commit/d953400d9c7f676fea70fee51264c331ba7e1067))
    - Python docs ([`10bf0df`](https://github.com/fiberseq/fibertools-rs/commit/10bf0df63acca5812e3d2d3ad5134bc75abe73b4))
    - Readme ([`cdde18e`](https://github.com/fiberseq/fibertools-rs/commit/cdde18e3b119ef36d57feed2f4740b7bde30c011))
    - I have a minimal working python package yay! ([`1f6fea3`](https://github.com/fiberseq/fibertools-rs/commit/1f6fea3cf07e1e798649df45d944a5396bf82902))
    - Start a python package ([`7e5fbde`](https://github.com/fiberseq/fibertools-rs/commit/7e5fbde9f0f52496fd1f2facbbd55d80d0b64a26))
    - Add to changelog ([`ac1987d`](https://github.com/fiberseq/fibertools-rs/commit/ac1987d33daf7fb7395ed1cf15b311390a15e60c))
    - Add to changelog ([`2f7c08e`](https://github.com/fiberseq/fibertools-rs/commit/2f7c08e2b24a07054d547088ac2c96463a6229d3))
    - Merge pull request #19 from fiberseq/development ([`202025d`](https://github.com/fiberseq/fibertools-rs/commit/202025dcdefdf88e0754ac78be03731c707c0b1b))
    - Update docs ([`c61762b`](https://github.com/fiberseq/fibertools-rs/commit/c61762b423129c7428ca946981febec491a9d853))
    - Update changelog ([`75a3a32`](https://github.com/fiberseq/fibertools-rs/commit/75a3a320729cfcedb3c99fd084cb2db16104ed70))
    - Unsed ([`4bc52a2`](https://github.com/fiberseq/fibertools-rs/commit/4bc52a229ef8be46eac45f765cfec6842d5d1d36))
    - Start a changelog ([`2aa57f5`](https://github.com/fiberseq/fibertools-rs/commit/2aa57f5fc1d7c38701964afea6ed5c20a76eb0d2))
    - Reorganize test data ([`966ac80`](https://github.com/fiberseq/fibertools-rs/commit/966ac80890e820da571164e01e360fc4c5507845))
    - Merge pull request #18 from fiberseq/development ([`2369f27`](https://github.com/fiberseq/fibertools-rs/commit/2369f2702c9caeb18d6c8475c9d1a52a6cbe423e))
</details>

## v0.2.5 (2023-07-18)

<csr-id-082190756233e077b692839f979df466fd5c7239/>

The default output for `ft extract` is now in reference coordinates instead of molecular.

`ft extract` can now produce compressed files by adding the gz extension. e.g.:

```
ft extract -t 16 --m6a m6a.bed.gz --cpg cpg.bed.gz --nuc nuc.bed.gz --msp msp.bed.gz --all all.bed.gz ../PS00243_ft.bam
```

And this is `bgzp` compression, so it is compatible with `tabix` indexing.

### Chore

- <csr-id-082190756233e077b692839f979df466fd5c7239/> Release fibertools-rs version 0.2.5

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 11 commits contributed to the release over the course of 6 calendar days.
 - 6 days passed between releases.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.2.5 ([`0821907`](https://github.com/fiberseq/fibertools-rs/commit/082190756233e077b692839f979df466fd5c7239))
    - Irnogre ([`938df0e`](https://github.com/fiberseq/fibertools-rs/commit/938df0e7e1aa88d05eb8b7e247eb88370e11d2b7))
    - Correct path ([`ac0166b`](https://github.com/fiberseq/fibertools-rs/commit/ac0166bb5cc1ac8b0b93df71ced2f8fd7689564f))
    - Better var name ([`fb3ad6d`](https://github.com/fiberseq/fibertools-rs/commit/fb3ad6d7c35e930cfb1ead67d494844eea62de93))
    - Fix off by one that sometimes happened in the final msp when doing nuc/msp calling ([`da5763d`](https://github.com/fiberseq/fibertools-rs/commit/da5763d24a322db2f1e6bdc53af57358776a51f4))
    - Fix off by one in minus strand nuc and msp coordiantes ([`91e2bc5`](https://github.com/fiberseq/fibertools-rs/commit/91e2bc558ab94b85e1ede25dccb2bf0e2878a023))
    - Merge pull request #17 from fiberseq/development ([`8f3d278`](https://github.com/fiberseq/fibertools-rs/commit/8f3d278d3285f8b5eafd429d5e6174b79adbbc49))
    - Check for bgz ([`14724dd`](https://github.com/fiberseq/fibertools-rs/commit/14724dd9270db0af9c12434fbc40252fd684b2fd))
    - Rewrite ft center so I can use threads? ([`7180f8e`](https://github.com/fiberseq/fibertools-rs/commit/7180f8e126a72fa49c4761e166bde1535c306aa2))
    - Rewrite ft center so I can use threads? ([`4485c5e`](https://github.com/fiberseq/fibertools-rs/commit/4485c5e46534b3bcdebf2295709e06830b404279))
    - Merge pull request #16 from fiberseq/development ([`759ee9f`](https://github.com/fiberseq/fibertools-rs/commit/759ee9fae9484f20255ae90b431587302c1f1cba))
</details>

## v0.2.4 (2023-07-11)

<csr-id-98a51ab3030ed03ce6adc7e160ee7feed49a7222/>

### Chore

- <csr-id-98a51ab3030ed03ce6adc7e160ee7feed49a7222/> Release fibertools-rs version 0.2.4

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 22 commits contributed to the release over the course of 13 calendar days.
 - 17 days passed between releases.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.2.4 ([`98a51ab`](https://github.com/fiberseq/fibertools-rs/commit/98a51ab3030ed03ce6adc7e160ee7feed49a7222))
    - Unsed deps ([`4603bb6`](https://github.com/fiberseq/fibertools-rs/commit/4603bb6920fd42e0c74d2b8d240d50b2d34c8a02))
    - Unsed deps ([`b542514`](https://github.com/fiberseq/fibertools-rs/commit/b542514c3353d5afba985763a51645a10cacdf1e))
    - Ci ([`60c065e`](https://github.com/fiberseq/fibertools-rs/commit/60c065e0515df493957929a0de16f916300df761))
    - Bio-io updates ([`8615452`](https://github.com/fiberseq/fibertools-rs/commit/8615452ba8743caa7b46c39f7acb583abe7ee129))
    - Merge pull request #15 from fiberseq/development ([`d8cfd67`](https://github.com/fiberseq/fibertools-rs/commit/d8cfd6773c37691c4aac73837a5124f9c7f4c6b7))
    - Change extract defaults ([`b009dda`](https://github.com/fiberseq/fibertools-rs/commit/b009ddacc3c07f9f86562ce25d987d04ab717c11))
    - CI ([`57880cc`](https://github.com/fiberseq/fibertools-rs/commit/57880cc16b29cbc250f9dfbf68dffd3c99f3c0fa))
    - Rename ([`a00e0d3`](https://github.com/fiberseq/fibertools-rs/commit/a00e0d30cf80d6ac987d4bfbcaa0bf46b128a067))
    - Set vars ([`eedb613`](https://github.com/fiberseq/fibertools-rs/commit/eedb613c2f4374ba826d1a8bca5402b97d5ccdee))
    - Fix version in bamwriter ([`c3e3781`](https://github.com/fiberseq/fibertools-rs/commit/c3e378118394b9d2deb8d7dde45ac7de28b48a3b))
    - Adding a compressed bgzip writer for files with .gz ([`e425ee3`](https://github.com/fiberseq/fibertools-rs/commit/e425ee32388b2870ae80e644af3f83d6ec167279))
    - Adding a buff reader that can be compressed optionally ([`4235b77`](https://github.com/fiberseq/fibertools-rs/commit/4235b7768e2530e0721cf70c9181811066ff9b75))
    - Making my own io package ([`d001c7f`](https://github.com/fiberseq/fibertools-rs/commit/d001c7f1e2cbd34bb5000093ba9565dc69ec3f65))
    - Help update ([`cfb85d9`](https://github.com/fiberseq/fibertools-rs/commit/cfb85d9a92af7765cf9921a474c6d1a1597abcd3))
    - Merge pull request #13 from Dylan-DPC/patch-1 ([`5bae1a1`](https://github.com/fiberseq/fibertools-rs/commit/5bae1a1195b14c7207a8cdff358fa23a3088f8f2))
    - Update Cargo.toml ([`14d0629`](https://github.com/fiberseq/fibertools-rs/commit/14d0629d03b98bbddc73f933158f472984aea34b))
    - Update docker to avoid double building in CI actions ([`19711d0`](https://github.com/fiberseq/fibertools-rs/commit/19711d0737e18354660c2a272d4b260a303a4a7f))
    - Merge pull request #12 from lvclark/patch-1 ([`dc8e230`](https://github.com/fiberseq/fibertools-rs/commit/dc8e2309953b0abcc51f5741be716202d2592a92))
    - Add HTSlib ([`fd058ba`](https://github.com/fiberseq/fibertools-rs/commit/fd058bac7084510f4c652c1da7e59df0db26bc77))
    - Put ft in a path that doesn't need root access ([`b024da2`](https://github.com/fiberseq/fibertools-rs/commit/b024da202fe250d89c59d4860639b9429ace8a36))
    - Uncomment compile line ([`d792d98`](https://github.com/fiberseq/fibertools-rs/commit/d792d98b46ec2b46c36ba3837a06db424c3acba7))
</details>

## v0.2.3 (2023-06-24)

<csr-id-6bb81e5a15f2f18bfe1a31b2ba207248e6a4e7b2/>

### Chore

- <csr-id-6bb81e5a15f2f18bfe1a31b2ba207248e6a4e7b2/> Release fibertools-rs version 0.2.3

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 20 commits contributed to the release over the course of 5 calendar days.
 - 8 days passed between releases.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.2.3 ([`6bb81e5`](https://github.com/fiberseq/fibertools-rs/commit/6bb81e5a15f2f18bfe1a31b2ba207248e6a4e7b2))
    - Update bamfile for release ([`5213be8`](https://github.com/fiberseq/fibertools-rs/commit/5213be8b2d0467c48086624cbd94b1a20573d2dc))
    - Fixing addition of nucleosomes to m6a in predict ([`7c50fab`](https://github.com/fiberseq/fibertools-rs/commit/7c50fabb8ccf2ae217d653fb642c23593964878e))
    - Add test ([`2545e26`](https://github.com/fiberseq/fibertools-rs/commit/2545e26b45521d19ecd13e0f01c9676a44cfc18a))
    - Add a command to srtip out basemods ([`8ecf75c`](https://github.com/fiberseq/fibertools-rs/commit/8ecf75cc9967d42c3125771e6c39704fcb77fef7))
    - Docker works :) ([`44016a7`](https://github.com/fiberseq/fibertools-rs/commit/44016a72c1dbbbff212ef206b7594319d053320d))
    - Update actions tool chain ([`24793fd`](https://github.com/fiberseq/fibertools-rs/commit/24793fd15a94bc3f64b074d91dfc62b98a568974))
    - Update actions tool chain ([`ee48a2d`](https://github.com/fiberseq/fibertools-rs/commit/ee48a2d266a1e7a95c8c8760383b17375e9b04af))
    - Update actions tool chain ([`3b15a91`](https://github.com/fiberseq/fibertools-rs/commit/3b15a9139bb898eb81aa858890865912fc8e0f89))
    - Ignore ([`2c1bcb4`](https://github.com/fiberseq/fibertools-rs/commit/2c1bcb46d432b2d5304adc4e909ef55ecab0ed98))
    - Fix docker env var ([`fc629e3`](https://github.com/fiberseq/fibertools-rs/commit/fc629e330af5dd620ba0d1ae85a66739bf3053d8))
    - Docker not needed here ([`4b89e52`](https://github.com/fiberseq/fibertools-rs/commit/4b89e52aaf5b264408f8cb21310afc685e26d30f))
    - Try to make gh actions use docker ([`57bb27d`](https://github.com/fiberseq/fibertools-rs/commit/57bb27d2bf8517aeaba9bce122a860d91151ab23))
    - Update docker ([`b28fe63`](https://github.com/fiberseq/fibertools-rs/commit/b28fe63de61c93f3a1c0d2785b460c0f161f9773))
    - Update docker ([`0e28b17`](https://github.com/fiberseq/fibertools-rs/commit/0e28b1739465f8ac791684d77b0125eec0b0c258))
    - Update docker ([`aba3f95`](https://github.com/fiberseq/fibertools-rs/commit/aba3f95b80caae822bff9ce97368243c97542690))
    - Update docker ([`1223538`](https://github.com/fiberseq/fibertools-rs/commit/1223538ac7c616950ebf6e6777548a4ddeca319d))
    - Simplify ([`5920f5b`](https://github.com/fiberseq/fibertools-rs/commit/5920f5bdc2af15812a5b335a5703be53dacc86ca))
    - Moving bamlift to its own crate ([`efd2a8b`](https://github.com/fiberseq/fibertools-rs/commit/efd2a8b6bf4e7ffe16ad3f91701cace194ea3967))
    - Moving bamlift to its own crate ([`ccd0969`](https://github.com/fiberseq/fibertools-rs/commit/ccd096963fe46b7ccbe00781b0da99cab87388ef))
</details>

## v0.2.2 (2023-06-15)

<csr-id-7cd3d115ffa9c824218ac8ea9f2454682cb91ca9/>

### Chore

- <csr-id-7cd3d115ffa9c824218ac8ea9f2454682cb91ca9/> Release fibertools-rs version 0.2.2

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 4 commits contributed to the release.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.2.2 ([`7cd3d11`](https://github.com/fiberseq/fibertools-rs/commit/7cd3d115ffa9c824218ac8ea9f2454682cb91ca9))
    - Add run tests and fix comments in bed files for center ([`8e00047`](https://github.com/fiberseq/fibertools-rs/commit/8e00047cdd769ad03a7e79197bfeb0e4bf00e6c4))
    - Add run tests ([`b5d1e8e`](https://github.com/fiberseq/fibertools-rs/commit/b5d1e8eca77d7fefec5eece5330f00a1e50abac9))
    - Hide some subcommands ([`03e7006`](https://github.com/fiberseq/fibertools-rs/commit/03e7006ab662e8fbb2ed0d3e11f7ffe1946ba4f1))
</details>

## v0.2.1 (2023-06-14)

<csr-id-53094a13be5fdc883bd87e46940bf93d48e41b1c/>

### Chore

- <csr-id-53094a13be5fdc883bd87e46940bf93d48e41b1c/> Release fibertools-rs version 0.2.1

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 5 commits contributed to the release over the course of 3 calendar days.
 - 3 days passed between releases.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.2.1 ([`53094a1`](https://github.com/fiberseq/fibertools-rs/commit/53094a13be5fdc883bd87e46940bf93d48e41b1c))
    - Readme update ([`24e26e8`](https://github.com/fiberseq/fibertools-rs/commit/24e26e857c3778c5e36c7a9afb0f68afba0b2d59))
    - Remove unused libs ([`fedd892`](https://github.com/fiberseq/fibertools-rs/commit/fedd892e1d360acb93b26f42a8b4d369c3fa8a48))
    - Bring back the colors ([`d336911`](https://github.com/fiberseq/fibertools-rs/commit/d3369115667fc7887eaf786de9f78a9f5504399f))
    - Bring back the colors ([`0af601a`](https://github.com/fiberseq/fibertools-rs/commit/0af601a8b6f8905be3ec79a0cc6db96bba940e61))
</details>

## v0.2.0 (2023-06-10)

<csr-id-cd387947c4aafd5d70d087650e1e5fe929808e81/>

### Chore

- <csr-id-cd387947c4aafd5d70d087650e1e5fe929808e81/> Release fibertools-rs version 0.2.0

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 19 commits contributed to the release over the course of 51 calendar days.
 - 52 days passed between releases.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.2.0 ([`cd38794`](https://github.com/fiberseq/fibertools-rs/commit/cd387947c4aafd5d70d087650e1e5fe929808e81))
    - Update version ([`d8e67d9`](https://github.com/fiberseq/fibertools-rs/commit/d8e67d9739ab0ad205050ffcef7ce2606854c80d))
    - Update version ([`6f541d3`](https://github.com/fiberseq/fibertools-rs/commit/6f541d31aad21e1b064203d10f6eeb8bb9d93b2d))
    - Update version ([`1931652`](https://github.com/fiberseq/fibertools-rs/commit/19316529bd0b427ab8da1007e7a6b263688ecd4e))
    - Update tch version ([`49983c3`](https://github.com/fiberseq/fibertools-rs/commit/49983c3ac085900397c505a702037790a28e9d50))
    - Add HP tag ([`3c845b4`](https://github.com/fiberseq/fibertools-rs/commit/3c845b40e17642c60c36774cbed8b1ddfcfc70cf))
    - Add HP tag ([`a6787ba`](https://github.com/fiberseq/fibertools-rs/commit/a6787ba6728d8dbc3d60b003ef4a9b34bf3ce954))
    - Add HP tag ([`9895f1c`](https://github.com/fiberseq/fibertools-rs/commit/9895f1c10c85fd32a9f12da1f8d59c2bae300447))
    - Add HP tag ([`4c78005`](https://github.com/fiberseq/fibertools-rs/commit/4c780051e44ee3f35a5539004f5ad03973638b48))
    - Add HP tag ([`6d9b74d`](https://github.com/fiberseq/fibertools-rs/commit/6d9b74d8c3d632fcde1e47242e33db00fb4a18d0))
    - Update for new revio model ([`08fcd9e`](https://github.com/fiberseq/fibertools-rs/commit/08fcd9e5d1890b0aa2c4243d9ddd097cf02adfae))
    - Cite ([`ce139ac`](https://github.com/fiberseq/fibertools-rs/commit/ce139ac5333ac1db84ab80d65b2f138301537d2e))
    - Fixing ignored min ml option ([`71f7799`](https://github.com/fiberseq/fibertools-rs/commit/71f77997a68dbb2b69f802375c49fba4ed4d865b))
    - Updating to new revio models ([`b0e7539`](https://github.com/fiberseq/fibertools-rs/commit/b0e753999c0885b0450c9d0bd4b42520380820da))
    - Updating to new revio models ([`b127d80`](https://github.com/fiberseq/fibertools-rs/commit/b127d8033e0d4a489a98009745d88bb791e31674))
    - Make logging more clear ([`daefc06`](https://github.com/fiberseq/fibertools-rs/commit/daefc067e979bfc95550149da557a3cb11013ef3))
    - Make logging more clear ([`ed8785a`](https://github.com/fiberseq/fibertools-rs/commit/ed8785ab203c7566426f220e2f6095e7bd6078b2))
    - Readme ([`2e25423`](https://github.com/fiberseq/fibertools-rs/commit/2e2542385974f2e584fcfe42b0bd76dedfefa6f3))
    - Readme ([`715e25d`](https://github.com/fiberseq/fibertools-rs/commit/715e25d9bbe82b36910984944b9ff3108a8f3bf3))
</details>

## v0.1.4 (2023-04-19)

<csr-id-4f95956d0a60a75b6b1a4d446829935426b2ffaf/>

### Chore

- <csr-id-4f95956d0a60a75b6b1a4d446829935426b2ffaf/> Release fibertools-rs version 0.1.4

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 44 commits contributed to the release over the course of 3 calendar days.
 - 36 days passed between releases.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.1.4 ([`4f95956`](https://github.com/fiberseq/fibertools-rs/commit/4f95956d0a60a75b6b1a4d446829935426b2ffaf))
    - Simplify ([`288ff21`](https://github.com/fiberseq/fibertools-rs/commit/288ff21d15508baaf035d962efad9702d7300190))
    - Simplify ([`687f41d`](https://github.com/fiberseq/fibertools-rs/commit/687f41d47ab003e6cd7cae5f7831d68fdc45e5a6))
    - Simplify ([`24f5930`](https://github.com/fiberseq/fibertools-rs/commit/24f5930a1e1827e3cdd7e1fb9533215c9ecf41b8))
    - New nucleosome caller is ready ([`2a856a7`](https://github.com/fiberseq/fibertools-rs/commit/2a856a796324bfe51934e9eaf70d945508a2f298))
    - New nucleosome caller is ready ([`6a5a370`](https://github.com/fiberseq/fibertools-rs/commit/6a5a3705457bce06f44ab149bed291f5a98c8add))
    - New nucleosome caller is ready ([`b4ee350`](https://github.com/fiberseq/fibertools-rs/commit/b4ee35093d0b0b8d7a2076feacd22fe074db8e2e))
    - New nucleosome caller is ready ([`39f045d`](https://github.com/fiberseq/fibertools-rs/commit/39f045da8c231a4ea4af30a1d03f71a408465575))
    - New nucleosome caller is ready ([`387fc4e`](https://github.com/fiberseq/fibertools-rs/commit/387fc4eab0613f7c47f9889789d26c44611c05b9))
    - New nucleosome caller is ready ([`ceb371c`](https://github.com/fiberseq/fibertools-rs/commit/ceb371c0a8c3401a0a4439f394bb427f71c13f15))
    - Skip double combine ([`b33dcae`](https://github.com/fiberseq/fibertools-rs/commit/b33dcaea2febcd99e5e4ce4d4d7f793380048f5e))
    - Skip double combine ([`80f421c`](https://github.com/fiberseq/fibertools-rs/commit/80f421c1ce3831bdf709752804506d3d4b6553be))
    - Skip double combine ([`330361c`](https://github.com/fiberseq/fibertools-rs/commit/330361c0fc96926730996ededf5b679da052347d))
    - Skip double combine ([`c63ba23`](https://github.com/fiberseq/fibertools-rs/commit/c63ba23b3bfa8a88e39971a5e6f69a2da7731c2d))
    - Skip double combine ([`7f2f23b`](https://github.com/fiberseq/fibertools-rs/commit/7f2f23bc4b4bc87b148247046b208b73f416e487))
    - Comdine over solitary m6a more agressiby ([`5cce88a`](https://github.com/fiberseq/fibertools-rs/commit/5cce88abed5f4e1fc82f7aa98c0be5f270e75f69))
    - Add rg to center ([`15d9ff1`](https://github.com/fiberseq/fibertools-rs/commit/15d9ff145a92970daf927dfa88d0b3e5515d52b6))
    - Fix forward rev ([`d953cec`](https://github.com/fiberseq/fibertools-rs/commit/d953cec6f129a9c1bdd6f3a7b2aac399433b3578))
    - Make dseg optional ([`feace70`](https://github.com/fiberseq/fibertools-rs/commit/feace70c6d331bfddc951979be6561eee41780eb))
    - Make dseg optional ([`d67754e`](https://github.com/fiberseq/fibertools-rs/commit/d67754ea94552e5e419c21438f959f33fb638d1e))
    - Make dseg optional ([`8e6a19e`](https://github.com/fiberseq/fibertools-rs/commit/8e6a19e1febe2ce0b068f98735a8b0a7d1f076c3))
    - Make dseg optional ([`2a3971b`](https://github.com/fiberseq/fibertools-rs/commit/2a3971b63b1b7091b12d15d4cd4b6ba769a91932))
    - Make dseg optional ([`a3291b9`](https://github.com/fiberseq/fibertools-rs/commit/a3291b9bfaef32cb9182ba57d8e0c0cacce44f27))
    - Make dseg optional ([`59f3f9f`](https://github.com/fiberseq/fibertools-rs/commit/59f3f9fde5ae717de0e295e4ec25da7cbae9b5f2))
    - Make dseg optional ([`ef5c7d2`](https://github.com/fiberseq/fibertools-rs/commit/ef5c7d2903809f04193705edc107859e7b8d8254))
    - Make dseg optional ([`e38363c`](https://github.com/fiberseq/fibertools-rs/commit/e38363c0d113ae31493550839f5414e301679aa9))
    - Make dseg optional ([`c5e42ed`](https://github.com/fiberseq/fibertools-rs/commit/c5e42ed3a6c0ffa14cb755bcf67aef3fa6adef5b))
    - D-seg paramenter on cmd line ([`eba6409`](https://github.com/fiberseq/fibertools-rs/commit/eba640965823d1c9650697195918a655971bb1a9))
    - Try dsegment ([`56cce35`](https://github.com/fiberseq/fibertools-rs/commit/56cce355483a8f92b790a23fd2c4ab76fbe5c9dc))
    - Update docs ([`404d60e`](https://github.com/fiberseq/fibertools-rs/commit/404d60edc217a49846a375a58916115fdde04bf4))
    - Improve help ([`b2e8f3e`](https://github.com/fiberseq/fibertools-rs/commit/b2e8f3e12aecf32253b164540a109726254ce70f))
    - Move nuc options to cli for predict-m6a ([`5c6c17b`](https://github.com/fiberseq/fibertools-rs/commit/5c6c17b502c679db52ae4d5d300375988911ca41))
    - Test update ([`4cfc00a`](https://github.com/fiberseq/fibertools-rs/commit/4cfc00aa14d1a1ba8c34b7b6178151d8286d878d))
    - Options change ([`ae86e48`](https://github.com/fiberseq/fibertools-rs/commit/ae86e4803d178d021ca1b59db0dccaf13e108142))
    - Add cli arguments ([`6530b50`](https://github.com/fiberseq/fibertools-rs/commit/6530b50522ced29ddb8cc0028df0177a67a3146e))
    - Test negatives ([`5d41679`](https://github.com/fiberseq/fibertools-rs/commit/5d41679c825e8a38cbdb3451e49c9e1c1af0a5f4))
    - Test negatives ([`86ea24c`](https://github.com/fiberseq/fibertools-rs/commit/86ea24cf3a899d30974bcf6c7b235238678fc357))
    - Test negatives ([`8acbfaf`](https://github.com/fiberseq/fibertools-rs/commit/8acbfaf0c31425b55ce889bbff8203c7f2622a6b))
    - Test negatives ([`b5cae9b`](https://github.com/fiberseq/fibertools-rs/commit/b5cae9b4b84fefd252a9be8e6ef81bc14df3aed2))
    - Test negatives ([`d17ccf2`](https://github.com/fiberseq/fibertools-rs/commit/d17ccf2896942f67b92dbd2dcbd160be64786583))
    - Test negatives ([`a178975`](https://github.com/fiberseq/fibertools-rs/commit/a178975974d8562cdc5d20974967b9d642ba26d6))
    - Adding a command that lets me add simple nuc calls ([`819a90b`](https://github.com/fiberseq/fibertools-rs/commit/819a90bd1255bbd8298ee60fde55b16a1bfccfe2))
    - Adding a command that lets me add simple nuc calls ([`33e0ee8`](https://github.com/fiberseq/fibertools-rs/commit/33e0ee84f8457d85ea52f45a430a61fcdc20dd8c))
    - Adding a command that lets me add simple nuc calls ([`cd3ec8f`](https://github.com/fiberseq/fibertools-rs/commit/cd3ec8fb718cf9353f43348b2b2c9741117f96d1))
</details>

## v0.1.3 (2023-03-14)

<csr-id-84e4ee9e1ae31396f8743d9d02c1af02dbe84bb7/>

### Chore

- <csr-id-84e4ee9e1ae31396f8743d9d02c1af02dbe84bb7/> Release fibertools-rs version 0.1.3

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 10 commits contributed to the release over the course of 16 calendar days.
 - 29 days passed between releases.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.1.3 ([`84e4ee9`](https://github.com/fiberseq/fibertools-rs/commit/84e4ee9e1ae31396f8743d9d02c1af02dbe84bb7))
    - Docs ([`2d97f26`](https://github.com/fiberseq/fibertools-rs/commit/2d97f26b66bd8f5e08fc66daacf74a925ff16677))
    - Ign ([`ae94728`](https://github.com/fiberseq/fibertools-rs/commit/ae94728e602bedfd1fdabb636078cd97e5aed968))
    - Adding 3.2 models ([`1193a76`](https://github.com/fiberseq/fibertools-rs/commit/1193a76e3ea8ac5f2225253d96ae770ebc0f7a07))
    - Further clean model chem seelection ([`32530e5`](https://github.com/fiberseq/fibertools-rs/commit/32530e50cbb3a4e23d1a304cc16b88c94c55b4a1))
    - Simplify model loading into the predict options ([`b60bf8b`](https://github.com/fiberseq/fibertools-rs/commit/b60bf8b8d43bb154a3f93be48f03d022e68238a3))
    - Try adding a reference option to ft center ([`ca015e0`](https://github.com/fiberseq/fibertools-rs/commit/ca015e01d5c6db82a982d919c3edbc60fc5af375))
    - Add the option to read the ML models from env variables FT_MODEL and FT_JSON ([`8314d8e`](https://github.com/fiberseq/fibertools-rs/commit/8314d8e62b21efe929c4e16c7509e7dd2bc312e4))
    - Add RG to extract ([`355088e`](https://github.com/fiberseq/fibertools-rs/commit/355088e59eac4c726b1d30212da3cc46d8bc59c2))
    - Theads help msg ([`595db94`](https://github.com/fiberseq/fibertools-rs/commit/595db94ef3b844ce49d829c25d79ee4b6d97f44d))
</details>

## v0.1.2 (2023-02-13)

<csr-id-7df4056ecc335f5be92f5688215a156623f927cd/>

### Chore

- <csr-id-7df4056ecc335f5be92f5688215a156623f927cd/> Release fibertools-rs version 0.1.2

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 3 commits contributed to the release over the course of 9 calendar days.
 - 10 days passed between releases.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.1.2 ([`7df4056`](https://github.com/fiberseq/fibertools-rs/commit/7df4056ecc335f5be92f5688215a156623f927cd))
    - Remvoe libtorch check since it is packaged in conda now ([`28c3614`](https://github.com/fiberseq/fibertools-rs/commit/28c3614ee3913451b3eb81697f904f104160930b))
    - Update url for new repo ([`021ae38`](https://github.com/fiberseq/fibertools-rs/commit/021ae389cf520391b5d3ccb066ffa8952ea32997))
</details>

## v0.1.1 (2023-02-02)

<csr-id-9f8f35b4f75ae3a944d137a8fd7bb562529a634e/>
<csr-id-d8f2ea597075a4f13be9d1d68c173c2476c02f63/>

### Chore

- <csr-id-9f8f35b4f75ae3a944d137a8fd7bb562529a634e/> Release fibertools-rs version 0.1.1
- <csr-id-d8f2ea597075a4f13be9d1d68c173c2476c02f63/> Release fibertools-rs version 0.1.1-alpha.4

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 6 commits contributed to the release.
 - 2 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.1.1 ([`9f8f35b`](https://github.com/fiberseq/fibertools-rs/commit/9f8f35b4f75ae3a944d137a8fd7bb562529a634e))
    - Bioconda now working ([`583ed05`](https://github.com/fiberseq/fibertools-rs/commit/583ed054e2b2cc5233f80e98e44c2fcabf8cec14))
    - Bioconda now working ([`d0f4943`](https://github.com/fiberseq/fibertools-rs/commit/d0f494361fa3f39a5d0cf85e1e0148ab27e4950a))
    - Fix ci ([`2d6fa56`](https://github.com/fiberseq/fibertools-rs/commit/2d6fa5617c814b651ee1f0e8f7214696c772d271))
    - Release fibertools-rs version 0.1.1-alpha.4 ([`d8f2ea5`](https://github.com/fiberseq/fibertools-rs/commit/d8f2ea597075a4f13be9d1d68c173c2476c02f63))
    - Restor version of ubunut ([`6ace556`](https://github.com/fiberseq/fibertools-rs/commit/6ace556673fe953e3bbdc50fbe9abff59e59aa0a))
</details>

## v0.1.1-alpha.3 (2023-02-01)

<csr-id-c1161138edd4b1bd7153a1cc5b724bb2eef29a4c/>

### Chore

- <csr-id-c1161138edd4b1bd7153a1cc5b724bb2eef29a4c/> Release fibertools-rs version 0.1.1-alpha.3

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 2 commits contributed to the release.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.1.1-alpha.3 ([`c116113`](https://github.com/fiberseq/fibertools-rs/commit/c1161138edd4b1bd7153a1cc5b724bb2eef29a4c))
    - Older version of ubunut ([`a043398`](https://github.com/fiberseq/fibertools-rs/commit/a04339858db5f2944f3a46be44eb2c3c367f1939))
</details>

## v0.1.1-alpha.2 (2023-02-01)

<csr-id-dda2ea4236a46e4cba2ac0f7415a871c09f5777a/>
<csr-id-c6f50df845bcd931b61cfa8151e4d38b1c9ccbdd/>

### Chore

- <csr-id-dda2ea4236a46e4cba2ac0f7415a871c09f5777a/> Release fibertools-rs version 0.1.1-alpha.2

### Chore

- <csr-id-c6f50df845bcd931b61cfa8151e4d38b1c9ccbdd/> Release fibertools-rs version 0.1.1-alpha.1

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 4 commits contributed to the release.
 - 2 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fibertools-rs version 0.1.1-alpha.2 ([`dda2ea4`](https://github.com/fiberseq/fibertools-rs/commit/dda2ea4236a46e4cba2ac0f7415a871c09f5777a))
    - Adding libtorch enviorment checks ([`3d0c77c`](https://github.com/fiberseq/fibertools-rs/commit/3d0c77c3bc89b1e6cce43e27ef1e92190339dcb2))
    - Release fibertools-rs version 0.1.1-alpha.1 ([`c6f50df`](https://github.com/fiberseq/fibertools-rs/commit/c6f50df845bcd931b61cfa8151e4d38b1c9ccbdd))
    - Change targets ([`1ad406a`](https://github.com/fiberseq/fibertools-rs/commit/1ad406a1de44d64ea84002442a9cab1d269d1ef8))
</details>

## v0.1.1-alpha.1 (2023-02-01)

<csr-id-c6f50df845bcd931b61cfa8151e4d38b1c9ccbdd/>

### Chore

- <csr-id-c6f50df845bcd931b61cfa8151e4d38b1c9ccbdd/> Release fibertools-rs version 0.1.1-alpha.1

## v0.1.1-alpha (2023-02-01)

<csr-id-3eede8cd716b8a03090a6f5dcc593460c6a943a1/>

### Chore

- <csr-id-3eede8cd716b8a03090a6f5dcc593460c6a943a1/> Release fibertools-rs version 0.1.1-alpha

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 6 commits contributed to the release.
 - 1 day passed between releases.
 - 1 commit was understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Change precommit ([`0b00d2d`](https://github.com/fiberseq/fibertools-rs/commit/0b00d2d69814fa25f5ecd0520e720b4b52fe81f0))
    - Release fibertools-rs version 0.1.1-alpha ([`3eede8c`](https://github.com/fiberseq/fibertools-rs/commit/3eede8cd716b8a03090a6f5dcc593460c6a943a1))
    - Default use cnn ([`9946f0f`](https://github.com/fiberseq/fibertools-rs/commit/9946f0fc7dc875392cd6ef55a537c4ead2f80970))
    - Testing a release manager ([`723bc6f`](https://github.com/fiberseq/fibertools-rs/commit/723bc6f0d7f7f7d808b041a6f6b8c50b6b2dd1db))
    - Testing a release manager ([`faf0620`](https://github.com/fiberseq/fibertools-rs/commit/faf062072eba594d006c2b4c50cea3b898193957))
    - Hide prediction under a new predict feature instead of partically under cnn ([`f86e1cf`](https://github.com/fiberseq/fibertools-rs/commit/f86e1cf628a3c3cb4db9a6ecd297467137d8dd4a))
</details>

## v0.1.0 (2023-01-30)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 134 commits contributed to the release over the course of 52 calendar days.
 - 52 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Make test smaller ([`45196f1`](https://github.com/fiberseq/fibertools-rs/commit/45196f1a62b9e7d0e007b6f5bd8729164dc729a9))
    - V0.1 ([`3f65061`](https://github.com/fiberseq/fibertools-rs/commit/3f65061287fd64e5176de2daff4466cdcc774ee9))
    - Check ([`14bb95f`](https://github.com/fiberseq/fibertools-rs/commit/14bb95f9bf4ca2f6a8369631369d22c3fe790dfb))
    - Fix default model and model setting flags in clap ([`a1bc54d`](https://github.com/fiberseq/fibertools-rs/commit/a1bc54d136647163f6ac7685b4093b4a28848338))
    - Change deafult to semi supervized cnn ([`ddfdc6e`](https://github.com/fiberseq/fibertools-rs/commit/ddfdc6e53f1620604d4197ec090bdd952e7229a3))
    - Fixing sequences reported by center that are on the reverse complement and check 6th column for strand and update docs ([`7bc77cc`](https://github.com/fiberseq/fibertools-rs/commit/7bc77ccdeeaef3269e24727ec735cb40cc2b0879))
    - Going back ([`e3a729e`](https://github.com/fiberseq/fibertools-rs/commit/e3a729ec15524fd4109cf98e50107fa852a60b37))
    - Add min ec coverage ([`1dbc5ed`](https://github.com/fiberseq/fibertools-rs/commit/1dbc5ed4badc14af6c7b476435acfe7b7bee6207))
    - Models ([`18f7634`](https://github.com/fiberseq/fibertools-rs/commit/18f7634b4ced91d78751f5b2eff95ababe6bcd63))
    - Models ([`7fe8cd4`](https://github.com/fiberseq/fibertools-rs/commit/7fe8cd45670fa31cf9c86b9d8e77775cf5bc032b))
    - Models ([`ceb21ef`](https://github.com/fiberseq/fibertools-rs/commit/ceb21efb7210d6f1340f3c5e0bd8381854a802ff))
    - Models ([`53ab799`](https://github.com/fiberseq/fibertools-rs/commit/53ab7994ec3602a09b9365ce69c4d5f8432673c2))
    - Try 225 initalization ([`909a5bc`](https://github.com/fiberseq/fibertools-rs/commit/909a5bc575667c16484f0faf375a416a011f1753))
    - Try 225 initalization ([`45f087b`](https://github.com/fiberseq/fibertools-rs/commit/45f087bff15f82b8b073c3045fd8f84c5526616d))
    - Try 225 initalization ([`3f7f202`](https://github.com/fiberseq/fibertools-rs/commit/3f7f202543c821cd93147170b04ac70e600ba819))
    - Try 225 initalization ([`61be875`](https://github.com/fiberseq/fibertools-rs/commit/61be8754f452b65ccd9b83cacb8d731512dc2354))
    - Try 225 initalization ([`d06f478`](https://github.com/fiberseq/fibertools-rs/commit/d06f4780bb04116e416c4910f6015c0c29bb5714))
    - Try 225 initalization ([`4d17b9d`](https://github.com/fiberseq/fibertools-rs/commit/4d17b9dd0037fdb36f56f8ee7ec0365da6dc3d10))
    - Try 225 initalization ([`2d59ffe`](https://github.com/fiberseq/fibertools-rs/commit/2d59ffecf0c42032316b954e994cca0b3a25d0c3))
    - Try 225 initalization ([`6837a01`](https://github.com/fiberseq/fibertools-rs/commit/6837a01aff9efdc3483b12f1cc89a253bd2537bb))
    - Try 225 initalization ([`b191f5f`](https://github.com/fiberseq/fibertools-rs/commit/b191f5f754dddf5c791692b334ac5765646882c0))
    - Try 225 initalization ([`f31840d`](https://github.com/fiberseq/fibertools-rs/commit/f31840d1b53deefde8f4c0e525f26b8251efbc1d))
    - Try 225 initalization ([`5c46adc`](https://github.com/fiberseq/fibertools-rs/commit/5c46adcfbb4978eb1b581d95027e947012934574))
    - Try 225 initalization ([`e66bdd0`](https://github.com/fiberseq/fibertools-rs/commit/e66bdd0989cd35ce4f8e9bc7d1a481de7d98a550))
    - Try 225 initalization ([`3837507`](https://github.com/fiberseq/fibertools-rs/commit/3837507d25124cebce728b630eb5e2ce7f295944))
    - Try 225 initalization ([`f920711`](https://github.com/fiberseq/fibertools-rs/commit/f920711d084c7e0a19257c6bde99fa5913cc16eb))
    - Try 225 initalization ([`2bcb25b`](https://github.com/fiberseq/fibertools-rs/commit/2bcb25bf43e95cc2136880655b19e9424bf77590))
    - Try 225 initalization ([`54fa3fb`](https://github.com/fiberseq/fibertools-rs/commit/54fa3fbea666be44f744374ed47abf3820ff2316))
    - Try 225 initalization ([`548655a`](https://github.com/fiberseq/fibertools-rs/commit/548655a897d04edbabe6563b305ff1d43e763629))
    - Try 225 initalization ([`e97ab5a`](https://github.com/fiberseq/fibertools-rs/commit/e97ab5a6aabe5ad282235b9bef32538a7c083328))
    - Try 225 initalization ([`729e9e4`](https://github.com/fiberseq/fibertools-rs/commit/729e9e4e23a94fb1f2f0450af3d8734033d3ed23))
    - Try 225 initalization ([`fbc82f5`](https://github.com/fiberseq/fibertools-rs/commit/fbc82f51d3c0499434b9febfaec2cb42885b7226))
    - Try 225 initalization ([`31926e6`](https://github.com/fiberseq/fibertools-rs/commit/31926e601c894a8605433bb2c23442bcca4c2d43))
    - Try 225 initalization ([`2086b51`](https://github.com/fiberseq/fibertools-rs/commit/2086b51c003b6ae30a34ebf55c169de40a7949ae))
    - Revio models ([`aebd08e`](https://github.com/fiberseq/fibertools-rs/commit/aebd08e94a6a0e53fd5acab43a9b0886aae2e672))
    - Revio models ([`37d99ec`](https://github.com/fiberseq/fibertools-rs/commit/37d99ec537ae8e0422ffb42aad8306fe456e9044))
    - Revio models ([`b891f68`](https://github.com/fiberseq/fibertools-rs/commit/b891f68ff9663c9efe550052aa554cd56d4e1c67))
    - Revio models ([`94a73fa`](https://github.com/fiberseq/fibertools-rs/commit/94a73fa1b46d24a98fe9b5bb6cc9206d172359cd))
    - Revio models ([`0cafca8`](https://github.com/fiberseq/fibertools-rs/commit/0cafca8076ccab29c7ac133e1750a6c911c95672))
    - Adding a subset of sequence to center ([`7a2f2f4`](https://github.com/fiberseq/fibertools-rs/commit/7a2f2f4e6d655206c34280ad9f8072bfe4aef7d6))
    - Adding batch of 10 reads all at once to gpu ([`c1947d5`](https://github.com/fiberseq/fibertools-rs/commit/c1947d5b2f6c36f9de6ef6189e9d586472c56efc))
    - Clean progress bar ([`9349850`](https://github.com/fiberseq/fibertools-rs/commit/9349850fcdb58a57b0c3014240daef23b02767e7))
    - Progress bar ([`f8763fb`](https://github.com/fiberseq/fibertools-rs/commit/f8763fb327695c32880662dfbda5cc346c9e3e07))
    - Add a clear kinetics command ([`c4d0103`](https://github.com/fiberseq/fibertools-rs/commit/c4d0103000a803784a85f5e358cceb25586827f0))
    - Add a optional distance filter to ft center ([`366ccc6`](https://github.com/fiberseq/fibertools-rs/commit/366ccc68ee06855b38dfa3870743d6c7f6503602))
    - Doc update ([`08b0878`](https://github.com/fiberseq/fibertools-rs/commit/08b08782c7ed714154c4ab339112503cdec48975))
    - Try next pytorch version ([`8d2a822`](https://github.com/fiberseq/fibertools-rs/commit/8d2a8225fc24834ff0a386749acc2494aeef718f))
    - Merge pull request #10 from mrvollger/dev ([`e3ef1e6`](https://github.com/fiberseq/fibertools-rs/commit/e3ef1e6411a0f603ea71f998493a8940e36dddd6))
    - Try next pytorch version ([`2e4eaca`](https://github.com/fiberseq/fibertools-rs/commit/2e4eaca457ad8ae6eb9a68ba94c07a8c7e7e8403))
    - Static link? ([`f68a4ec`](https://github.com/fiberseq/fibertools-rs/commit/f68a4ec59a95bc3bd263fd4843ac578b99df790c))
    - Adding a ml filter ([`839a4d2`](https://github.com/fiberseq/fibertools-rs/commit/839a4d24287357d1fe80c87e251212cf50f7c631))
    - Adding pb u16 converter ([`5f8d877`](https://github.com/fiberseq/fibertools-rs/commit/5f8d877700d3ec4fcdaf9905c4c262d0b5398c22))
    - Final (fingers crossed) values for ml filters for 2.0 and 2.2 chem ([`ee55a17`](https://github.com/fiberseq/fibertools-rs/commit/ee55a1711013061bf41d36d9824592ceba1266f4))
    - Update tables for models ([`265dce1`](https://github.com/fiberseq/fibertools-rs/commit/265dce1ddc59c40eea80c5c332b41a8ef62cdf0d))
    - Models out of date wrt tables ([`3d25480`](https://github.com/fiberseq/fibertools-rs/commit/3d25480d1793a9d13df1baa14e3d8258f777a2c1))
    - Read in a json prediction conversion table and apply it to the semi suppervized approach ([`8a57047`](https://github.com/fiberseq/fibertools-rs/commit/8a5704715f7deb291ac5799212f0ca81adc88af7))
    - Fmt ([`2e93877`](https://github.com/fiberseq/fibertools-rs/commit/2e93877a1fcbf95113e1e7f60e74029f8d30c7d4))
    - Fmt ([`fb78c73`](https://github.com/fiberseq/fibertools-rs/commit/fb78c730918e870eed40d9056e3b1fe9d0423c42))
    - Framework for converting floats to u8 precision calls ([`95185fa`](https://github.com/fiberseq/fibertools-rs/commit/95185faa33f27d807386203a1ddea220c757c6fa))
    - Limit width ([`49b90a9`](https://github.com/fiberseq/fibertools-rs/commit/49b90a95bb3aefbe2b9f1c373e7415347281bcbd))
    - Logos ([`d3d7f23`](https://github.com/fiberseq/fibertools-rs/commit/d3d7f23328610598dee9fba64342391d6a1c5222))
    - Updating default values for models based on chemsitry ([`f706aef`](https://github.com/fiberseq/fibertools-rs/commit/f706aef82e67d32774a9d19899740b3327199864))
    - Fix overwirrten model and add new 2.2 model ([`d6b479d`](https://github.com/fiberseq/fibertools-rs/commit/d6b479dff23d054a8d1f88551bd91e4b98d8b0ef))
    - Fix overwirrten model and add new 2.2 model ([`efa2b51`](https://github.com/fiberseq/fibertools-rs/commit/efa2b51302482b11887e988e5f86e711daa96d9b))
    - Fix overwirrten model and add new 2.2 model ([`d8a1ce3`](https://github.com/fiberseq/fibertools-rs/commit/d8a1ce389824a757856eda4c03510c0c0cec3715))
    - Adjust suggested value for semi ([`95af6b9`](https://github.com/fiberseq/fibertools-rs/commit/95af6b9abdcedf5951f10728ed54db0d927705b7))
    - Update model ([`eca5448`](https://github.com/fiberseq/fibertools-rs/commit/eca5448da750ed163fb3cd082b31958a4c0b38f3))
    - Update model ([`660a1ac`](https://github.com/fiberseq/fibertools-rs/commit/660a1ac6fa18c0c619e2d9eb0c57e2f394ac99dc))
    - Update model ([`08b35bd`](https://github.com/fiberseq/fibertools-rs/commit/08b35bdfe5c06b0191c73c74036ddf421cd9b50a))
    - Update model ([`c58c252`](https://github.com/fiberseq/fibertools-rs/commit/c58c252845dcb185c8e204f720384caf9ef83ef3))
    - Update model ([`33df823`](https://github.com/fiberseq/fibertools-rs/commit/33df823c8588d0fd149bbfd8f7bd3dd6928b1cce))
    - Update model ([`fb0ef59`](https://github.com/fiberseq/fibertools-rs/commit/fb0ef598a8529c355c0af6e2a65fbc91931ce826))
    - Update model ([`42d7be0`](https://github.com/fiberseq/fibertools-rs/commit/42d7be0a5d5910a644686b0ae44f178a87eb800e))
    - Update model ([`0e709cc`](https://github.com/fiberseq/fibertools-rs/commit/0e709cc0d480788bb8fa6ab484f111a052c7ad3b))
    - Update model ([`7602cb6`](https://github.com/fiberseq/fibertools-rs/commit/7602cb6f03943d32084644981ea04d6531776657))
    - Update model ([`8642198`](https://github.com/fiberseq/fibertools-rs/commit/8642198dd0a9d0871db02accd6c0d3602645e99e))
    - Update model ([`263ff71`](https://github.com/fiberseq/fibertools-rs/commit/263ff71a8473a6fcd9b246f840ce347a4ea6c9fe))
    - Update model ([`2e84b17`](https://github.com/fiberseq/fibertools-rs/commit/2e84b174f8fe2671ca8e2cafd63f0013e51eb048))
    - Update model ([`c2c637e`](https://github.com/fiberseq/fibertools-rs/commit/c2c637e9adedaa31de28998e5db684365f8a20d8))
    - Update model ([`d86c1c9`](https://github.com/fiberseq/fibertools-rs/commit/d86c1c92a23148d92a27f4ffac3534005d8d5fb8))
    - Update model ([`1d5568e`](https://github.com/fiberseq/fibertools-rs/commit/1d5568e9ee7c55ae5a76a6b4d66fa4f031c2314d))
    - Update model ([`652992c`](https://github.com/fiberseq/fibertools-rs/commit/652992c2e200f5c0d37bfdc614e0625a3dfba431))
    - Update model ([`7c83e6f`](https://github.com/fiberseq/fibertools-rs/commit/7c83e6fc3d5b1265ede574f9cadfed43f9f92032))
    - Update model ([`d1a480a`](https://github.com/fiberseq/fibertools-rs/commit/d1a480a01dc9de5482cb7e30c85f1b20800f59f1))
    - Update model ([`9d57185`](https://github.com/fiberseq/fibertools-rs/commit/9d571856b63e49dcb1baa6ac950cba0612de5de1))
    - Update model ([`1532b69`](https://github.com/fiberseq/fibertools-rs/commit/1532b69e85ae81d60812b33153d4172af23e48c3))
    - Update model ([`aed4e39`](https://github.com/fiberseq/fibertools-rs/commit/aed4e3951e0c2f974088b67de4a290f23841cefa))
    - Update model ([`ace39e1`](https://github.com/fiberseq/fibertools-rs/commit/ace39e1799efb15821e0d11a022bd99ee16730c2))
    - Update model ([`3838aad`](https://github.com/fiberseq/fibertools-rs/commit/3838aad6941ad934d2a12156d8ff72ea31ca40fa))
    - Update model ([`1bed0db`](https://github.com/fiberseq/fibertools-rs/commit/1bed0dbf9b33ffe439571386960bd35680bd765d))
    - Update model ([`ab4a2af`](https://github.com/fiberseq/fibertools-rs/commit/ab4a2af0964a89b3d25ee5539f26d3dd834654fd))
    - Update model ([`6a093dc`](https://github.com/fiberseq/fibertools-rs/commit/6a093dc0a6f795f69492e2afb1747d123d4e1a97))
    - Update model ([`b137ca1`](https://github.com/fiberseq/fibertools-rs/commit/b137ca120c123d72dd8a21b5522da2431a17493a))
    - Update model ([`7123d61`](https://github.com/fiberseq/fibertools-rs/commit/7123d6170ddf39797d3860acfd3df11666743227))
    - Update model ([`83c0b2e`](https://github.com/fiberseq/fibertools-rs/commit/83c0b2ec82f601a0094e3a2ebb5e71f2aac56303))
    - Update model ([`383d09f`](https://github.com/fiberseq/fibertools-rs/commit/383d09f06e31dbbc18bad0c715eae04cfdd6e39f))
    - Update model ([`a801f7a`](https://github.com/fiberseq/fibertools-rs/commit/a801f7ab147edfb96c1d9f99a300cfa0a770d52e))
    - Update model ([`f6c1513`](https://github.com/fiberseq/fibertools-rs/commit/f6c15134b20d82fdf4854115c86d759aa2a7ec24))
    - No semi score ([`d2f58b2`](https://github.com/fiberseq/fibertools-rs/commit/d2f58b235c783989021783d30f946bdc0bfc28fe))
    - Custom semi adjust ([`0a9ef81`](https://github.com/fiberseq/fibertools-rs/commit/0a9ef81fc31c172a8c11b18628f39302a46a5955))
    - Custom semi adjust ([`0758311`](https://github.com/fiberseq/fibertools-rs/commit/07583118d82194caeb89914724156c4e3c9689b2))
    - Custom semi adjust ([`7bb553e`](https://github.com/fiberseq/fibertools-rs/commit/7bb553e76e17c56b9f29ae9abc9e6d5e676c219c))
    - Custom semi adjust ([`c0b3b73`](https://github.com/fiberseq/fibertools-rs/commit/c0b3b7303f45116ab3994ff1648ca9ce718156bf))
    - Custom semi adjust ([`25b8811`](https://github.com/fiberseq/fibertools-rs/commit/25b88110e6848b46cfbb5c134ca1d16a6f4da687))
    - Custom semi adjust ([`b28f2dc`](https://github.com/fiberseq/fibertools-rs/commit/b28f2dcedf5b3b550edc3caac5c5f21a558b135a))
    - Custom semi adjust ([`f7ddc51`](https://github.com/fiberseq/fibertools-rs/commit/f7ddc51372ea7ea5b9f0e5be3ec33e3ba73b38bc))
    - Custom semi adjust ([`29b9f2c`](https://github.com/fiberseq/fibertools-rs/commit/29b9f2c3163ac7be3523593286fa7408f09d85fd))
    - Custom semi adjust ([`de586ff`](https://github.com/fiberseq/fibertools-rs/commit/de586ff65d324a9f5c7c74eaa3e9997fa2a0a766))
    - Custom semi adjust ([`66fec59`](https://github.com/fiberseq/fibertools-rs/commit/66fec59a0e718566d90578fde7d3d9e52f86b939))
    - Custom semi adjust ([`db55c50`](https://github.com/fiberseq/fibertools-rs/commit/db55c5056d421f4006399b6356b6fd68ca28c1b5))
    - Custom semi adjust ([`27bb5d1`](https://github.com/fiberseq/fibertools-rs/commit/27bb5d1f413ca47179d984976800f79342490e80))
    - Custom semi adjust ([`758be9b`](https://github.com/fiberseq/fibertools-rs/commit/758be9b529643fc26248851284bf769756e942cf))
    - Custom semi adjust ([`a595d70`](https://github.com/fiberseq/fibertools-rs/commit/a595d70b8f4b7e8e1a2dddb0637aef1d28585729))
    - Custom semi adjust ([`6a1e6f1`](https://github.com/fiberseq/fibertools-rs/commit/6a1e6f1b106013b921cfe50333e47bb0d35370d6))
    - Larger full float ([`1a2f942`](https://github.com/fiberseq/fibertools-rs/commit/1a2f9425888140427a5555daa6d9e33030032866))
    - Restore full float if needed ([`5ec3156`](https://github.com/fiberseq/fibertools-rs/commit/5ec31564fdafaf5a1d9da956dda06ee869c4952c))
    - Sam flag docs and man page generation ([`29d21d9`](https://github.com/fiberseq/fibertools-rs/commit/29d21d9fe758df7e0d522889064896f1136b1c07))
    - Add completions as a subcommand ([`0cbc39a`](https://github.com/fiberseq/fibertools-rs/commit/0cbc39a48aea83a3f4c98733e74ed4dcb00b37e1))
    - Adding a sam flag, fixes #9 ([`6643419`](https://github.com/fiberseq/fibertools-rs/commit/66434192371259e42b872f0c1dd8647c1d38b20d))
    - Update 2.0 semi model ([`f9f1275`](https://github.com/fiberseq/fibertools-rs/commit/f9f127595becd27c9e26bccf9831e5a3632ec7a6))
    - Fix for unaligned ([`01fd136`](https://github.com/fiberseq/fibertools-rs/commit/01fd13615116a652712f25e2fc10da6f92b3b115))
    - Move target name into fiberdata so I can use a par iter, since HeaderView does not have sync... ([`a57548c`](https://github.com/fiberseq/fibertools-rs/commit/a57548cc03c10fcc8751176997364ecebe23f3d2))
    - Use bam iterator in ft extract ([`13b0b35`](https://github.com/fiberseq/fibertools-rs/commit/13b0b35f7f40e67c87675687ab14cd1f72118092))
    - Dont mess with gpu, nvm ([`574c66d`](https://github.com/fiberseq/fibertools-rs/commit/574c66dcbb9e9eda80501e76027c34b54a7fe96c))
    - Dont mess with gpu ([`f6dfd7f`](https://github.com/fiberseq/fibertools-rs/commit/f6dfd7fb710c905bf853c70d705db214356344d1))
    - Comment ([`46ceed8`](https://github.com/fiberseq/fibertools-rs/commit/46ceed8fbf39c94a515f7f98ccf54fb29cf5ff5c))
    - Make sure tch is only useing 1 thread so rayon can use all the other ones ([`bc45715`](https://github.com/fiberseq/fibertools-rs/commit/bc45715126e80688d92a6eb769dd5e68b3dbebb4))
    - Simplify predict chunk for cnn, maybe ([`0db8233`](https://github.com/fiberseq/fibertools-rs/commit/0db8233c134f78747a1192c609a52bca71613292))
    - Actually that slowed it down ([`8e4bfdf`](https://github.com/fiberseq/fibertools-rs/commit/8e4bfdf8a8366b2518167c79c3648906ff4ef860))
    - Speed up cpu cnn model by reducing batch size ([`8e1f47a`](https://github.com/fiberseq/fibertools-rs/commit/8e1f47a4b46cfea91eaebbcb8eeb6da3f4302b2a))
    - Fake model quickly ([`c24a678`](https://github.com/fiberseq/fibertools-rs/commit/c24a67862e1bcfd616e060d7ae81ef76d84a3743))
    - Inital adding of a m6A ML filter during prediction ([`95e17f0`](https://github.com/fiberseq/fibertools-rs/commit/95e17f0ccb71ec6a4d3ca819aeaddb31824d9b7f))
    - Inital adding of a m6A ML filter during prediction ([`d9dbcad`](https://github.com/fiberseq/fibertools-rs/commit/d9dbcad107291b8bf903b5d7908a226bfd0c638f))
    - Merge pull request #8 from mrvollger/main ([`6e28135`](https://github.com/fiberseq/fibertools-rs/commit/6e281357ed98c204f58b4f90d91960f06f3368df))
</details>

## v0.0.11 (2022-12-09)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 63 commits contributed to the release over the course of 17 calendar days.
 - 17 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Ignore ([`a427005`](https://github.com/fiberseq/fibertools-rs/commit/a42700559815dbee6ef047625065f41d11cb33d8))
    - Merge pull request #7 from mrvollger/dev ([`794a668`](https://github.com/fiberseq/fibertools-rs/commit/794a66898105fbaffc3b91a9a7a5086800c06ecb))
    - Merge pull request #6 from mrvollger/main ([`168d985`](https://github.com/fiberseq/fibertools-rs/commit/168d985e8c317080c76bc8249962d571278c5f14))
    - Doc update ([`5b85aff`](https://github.com/fiberseq/fibertools-rs/commit/5b85aff92564fd37e7074fefb98397b34dd5caa0))
    - Merge pull request #5 from mrvollger/dev ([`64f5111`](https://github.com/fiberseq/fibertools-rs/commit/64f51111c6e8c9a6be63fb9132124f8da682622c))
    - Change to install from main ([`4e8c53d`](https://github.com/fiberseq/fibertools-rs/commit/4e8c53de5fcbe0f53074d2cc5c8e65fab7d9eb66))
    - Nuc readme details ([`34395e7`](https://github.com/fiberseq/fibertools-rs/commit/34395e7e3e3f79c3376db163533a03ab59ddead8))
    - 2.2 model final ([`83e4d0d`](https://github.com/fiberseq/fibertools-rs/commit/83e4d0dba3e690ba63fe3f44e64601c4c7f85f07))
    - ML filters ([`7b1e479`](https://github.com/fiberseq/fibertools-rs/commit/7b1e479a81f800dac631c14d8cfcf4d99efb8c8f))
    - ML filters ([`02742c0`](https://github.com/fiberseq/fibertools-rs/commit/02742c0c619d3f47d46863612a8a345b92ed6a0d))
    - Todo ([`cc3c1ea`](https://github.com/fiberseq/fibertools-rs/commit/cc3c1eab598b4c5e4f03a98347a068d39c9e687a))
    - Semi model update ([`d4f8146`](https://github.com/fiberseq/fibertools-rs/commit/d4f8146e94940eaba350cc93e06c568fd38b8656))
    - Adding use of namedtempfile ([`663f296`](https://github.com/fiberseq/fibertools-rs/commit/663f2961441dde81f8a3e0e44c7166b3315b0146))
    - Update model ([`d3c2590`](https://github.com/fiberseq/fibertools-rs/commit/d3c259054dafe09b11bf4bedd3327ddab9bc9c48))
    - Trying new website themes ([`1ead9fd`](https://github.com/fiberseq/fibertools-rs/commit/1ead9fdfe807fd1beedc6cdcb1ebba9a3aa1e579))
    - Add logo ([`7fd71b2`](https://github.com/fiberseq/fibertools-rs/commit/7fd71b2531648106c6d2a8b938302d695d34b901))
    - Move images ([`005ea70`](https://github.com/fiberseq/fibertools-rs/commit/005ea70a53a8a83d3f49e5298eeb5b76b3e4f510))
    - Adding 2.2 model ([`00a4461`](https://github.com/fiberseq/fibertools-rs/commit/00a44610df1362f5635aa276e830ada54c4a2eeb))
    - Adding semi model file ([`1ffae74`](https://github.com/fiberseq/fibertools-rs/commit/1ffae74badbf89e67d5d320fc734f048e66e98be))
    - Logging ([`dad490d`](https://github.com/fiberseq/fibertools-rs/commit/dad490d3a78295dae2959c684450dfb8eb8b2256))
    - Add and change env vars ([`ee2a687`](https://github.com/fiberseq/fibertools-rs/commit/ee2a687543e47caf657fa22374516e0be0993ce5))
    - Fmt ([`4392f1f`](https://github.com/fiberseq/fibertools-rs/commit/4392f1f4b7eb5c005afc59aeb294a690db44f198))
    - Add version to git hash ([`3b19223`](https://github.com/fiberseq/fibertools-rs/commit/3b192236793fe62b20b4693876dcc296a0749821))
    - Add git hash for long version ([`54cb7ae`](https://github.com/fiberseq/fibertools-rs/commit/54cb7ae89bbc6e106807b96aa7df20c52766ff39))
    - Add option for semi supervized ([`a80df7b`](https://github.com/fiberseq/fibertools-rs/commit/a80df7ba665207ac9fe256733728d19054faf367))
    - Add colors ([`aef5a7e`](https://github.com/fiberseq/fibertools-rs/commit/aef5a7e77276e09f68500950841daf541f510b15))
    - Ci ([`297ff99`](https://github.com/fiberseq/fibertools-rs/commit/297ff9981b6a99e362371fcdacd93d30e847702c))
    - Docs ([`98dd885`](https://github.com/fiberseq/fibertools-rs/commit/98dd885d0bae12eda518968e8c301a957a5eefae))
    - Moving base mods and adding function to get forward bases ([`eac73c6`](https://github.com/fiberseq/fibertools-rs/commit/eac73c65bdacccb8aa612bd5cf753da5dacde327))
    - Description ([`76aba78`](https://github.com/fiberseq/fibertools-rs/commit/76aba78e1eae1591d2c526112adb01711b5584e5))
    - Ideas for v0.0.11 ([`1a9114d`](https://github.com/fiberseq/fibertools-rs/commit/1a9114d82fbd5dd1b65507c27b9e60976110c244))
    - Files moved ([`ff85512`](https://github.com/fiberseq/fibertools-rs/commit/ff855127f77673e445ac14b2fddcf4080e1c6031))
    - Docs ([`d3b8e62`](https://github.com/fiberseq/fibertools-rs/commit/d3b8e629449dd40764f2e2ebb58d89fd797f3cd9))
    - Merge pull request #4 from mrvollger/gpu ([`778ad7a`](https://github.com/fiberseq/fibertools-rs/commit/778ad7aafa0ad4d69808594f85b75c2fe6c57179))
    - Remove unused and check secondary ([`828e684`](https://github.com/fiberseq/fibertools-rs/commit/828e6840ddfb8483cfa233a3e4f25b9db9b33df7))
    - Fix sucsess counter ([`3688365`](https://github.com/fiberseq/fibertools-rs/commit/3688365655c95e68cc43109b6ff8bafb91fbee7c))
    - Move model shape into consts ([`37d8167`](https://github.com/fiberseq/fibertools-rs/commit/37d81677cc77c5ddb628089f86df0a05648b27d5))
    - Progress ([`db2b2e0`](https://github.com/fiberseq/fibertools-rs/commit/db2b2e0827559641342e5b8b6d81281b70994125))
    - Batchsize ([`7e2cd39`](https://github.com/fiberseq/fibertools-rs/commit/7e2cd39b5c1929b8a96c13595a47e4e5ddd5375f))
    - Msg ([`a4e96b0`](https://github.com/fiberseq/fibertools-rs/commit/a4e96b026ae9469468aa55da24a70dd68435fb9d))
    - Msg ([`ead03fa`](https://github.com/fiberseq/fibertools-rs/commit/ead03faf56aeee5ac8f49be40d6a283cea930ff1))
    - Msg ([`9c88c03`](https://github.com/fiberseq/fibertools-rs/commit/9c88c03690a694fe19c5d39b9218a4ac76960e41))
    - Change xgb batch size to 1 because it is fastest on cpu ([`5f83ef4`](https://github.com/fiberseq/fibertools-rs/commit/5f83ef4a609c57cf43abbc6655cc66c819d4f355))
    - Add zero check ([`3a8c74b`](https://github.com/fiberseq/fibertools-rs/commit/3a8c74b74fd3c7073149042a5be9f9dd02661006))
    - Make batch size defaults ([`c45eeda`](https://github.com/fiberseq/fibertools-rs/commit/c45eeda30cdf9eb38fefbd51ffc71d61332cf24a))
    - Make batch size an option ([`5b5c86d`](https://github.com/fiberseq/fibertools-rs/commit/5b5c86d7c3298bbe931373389cb946aaed07a500))
    - Chunk ([`bb79e62`](https://github.com/fiberseq/fibertools-rs/commit/bb79e62b4e069aab8d4240b2702b16deecf40cd2))
    - Move progress ([`5d09a54`](https://github.com/fiberseq/fibertools-rs/commit/5d09a54d5d2d1f47e3b3cda8b764fa7974867db8))
    - Adding batch of 10 reads all at once to gpu ([`34c8cbe`](https://github.com/fiberseq/fibertools-rs/commit/34c8cbe6ae75e9c46e390849075cd33ca4940729))
    - Adding more install instructions for gpu ([`66d0eae`](https://github.com/fiberseq/fibertools-rs/commit/66d0eae9a7a86e6d05180506e88f305dac575abb))
    - Comment ([`ec4c1e9`](https://github.com/fiberseq/fibertools-rs/commit/ec4c1e990d87fa63f1bad1d0d07ae9f085a0ec4f))
    - Merge pull request #3 from mrvollger/gpu ([`1529b9e`](https://github.com/fiberseq/fibertools-rs/commit/1529b9ea32c01004a8e458bac4eb6a74c5194cf0))
    - Move prgress to after map ([`765ba52`](https://github.com/fiberseq/fibertools-rs/commit/765ba52016a8d43d2b508e091e7b944e08dff4c1))
    - Move model to gpu? ([`a928875`](https://github.com/fiberseq/fibertools-rs/commit/a928875bb494e27a65a1b6d178f70923fc9e446f))
    - Move vec to gpu? ([`b5853d5`](https://github.com/fiberseq/fibertools-rs/commit/b5853d5f45ea67d144476e78e05488e3e23496a1))
    - Move vec to gpu? ([`a213c7c`](https://github.com/fiberseq/fibertools-rs/commit/a213c7c03d5163a701246ce5dbc630db4337ebf9))
    - Check device ([`77d5cac`](https://github.com/fiberseq/fibertools-rs/commit/77d5cac9bed0537d8045867b29de13cbbc56aa78))
    - Remove dash ([`31728e6`](https://github.com/fiberseq/fibertools-rs/commit/31728e688be9831c62b770f377c8be7bf5fc1874))
    - Reorder readme ([`79c4974`](https://github.com/fiberseq/fibertools-rs/commit/79c497452e3c52a7ffef0c781ce23a1ad67a5831))
    - More details on cnn install ([`5bcbef8`](https://github.com/fiberseq/fibertools-rs/commit/5bcbef8ee5d9db7626dafcf82b4782053f67c8e7))
    - Fix program line for ft ([`c30ffa3`](https://github.com/fiberseq/fibertools-rs/commit/c30ffa3bdcbada6d2932656bae72e82dd033cf51))
    - Caps ([`ac9ecfb`](https://github.com/fiberseq/fibertools-rs/commit/ac9ecfbbb4ed22d423954348419ae48f50c98c88))
    - Clean ([`93839d3`](https://github.com/fiberseq/fibertools-rs/commit/93839d3e6667966fb51a9a9e7c3924d9ce84b6f5))
</details>

## v0.0.10 (2022-11-22)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 9 commits contributed to the release.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Clean ([`3f5fbec`](https://github.com/fiberseq/fibertools-rs/commit/3f5fbec24ce2e63ec7517d5abafb011a245f1072))
    - Clean ([`25d11b0`](https://github.com/fiberseq/fibertools-rs/commit/25d11b025c3e907c5f0d11a3e2a28628bd4742f6))
    - Get ready for the next version ([`540b76c`](https://github.com/fiberseq/fibertools-rs/commit/540b76c1e8789f7323fb46f42637ec3977dcb2bb))
    - Typo ([`0fc22d1`](https://github.com/fiberseq/fibertools-rs/commit/0fc22d1d6f5fb7715f0d5adb48515815eadd24fc))
    - Add doc strings ([`bcdd7d0`](https://github.com/fiberseq/fibertools-rs/commit/bcdd7d0ec4a077693fe117ecb7eb3a4dd73cb45a))
    - Docs ([`ddd998f`](https://github.com/fiberseq/fibertools-rs/commit/ddd998fe135da32f5e6ffa7b85d5f8a2e34c2de8))
    - Adding a PG tag to the m6a output ([`3d004aa`](https://github.com/fiberseq/fibertools-rs/commit/3d004aaffbac11b156e2ea67ea8a0abb99e6a770))
    - Typo ([`14d69aa`](https://github.com/fiberseq/fibertools-rs/commit/14d69aab04fc180cf158e8b0d30b7f70a64dc6d7))
    - Update release script for CNN ([`5c25138`](https://github.com/fiberseq/fibertools-rs/commit/5c25138922fd06b11e4f91509d2b336242cabe35))
</details>

## v0.0.9 (2022-11-21)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 15 commits contributed to the release.
 - 1 day passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - A new release for bioconda which will hopefully allow the cnn by default ([`85305ab`](https://github.com/fiberseq/fibertools-rs/commit/85305ab668c9298008b7d8c7116c88f108b02218))
    - Increase chunk size ([`a028151`](https://github.com/fiberseq/fibertools-rs/commit/a028151da378709cc81085417a145df2e53d7982))
    - For #2 try 2 ([`950d7db`](https://github.com/fiberseq/fibertools-rs/commit/950d7db10a84b9c55f4c53f854ac84c91e8b06bb))
    - For #2 ([`de08a90`](https://github.com/fiberseq/fibertools-rs/commit/de08a90f752b1a94b15486bdced4a42b6a093d21))
    - Stadardize message formating ([`178ef3c`](https://github.com/fiberseq/fibertools-rs/commit/178ef3cadd5922fa65b56398bd41b5700745d25e))
    - Stadardize message formating ([`0fc4761`](https://github.com/fiberseq/fibertools-rs/commit/0fc4761a6892e862e3dbb1c5c0e96677e1250f4e))
    - Add conda-forge to channel ([`dd65ce9`](https://github.com/fiberseq/fibertools-rs/commit/dd65ce9fab9ddc5e9ff9035eb60fd1653508b73d))
    - Add detail to cli ([`e081131`](https://github.com/fiberseq/fibertools-rs/commit/e081131fe8643caf0b2050f19c31a1ef58517f07))
    - Adding more headings in CLI ([`c8922fd`](https://github.com/fiberseq/fibertools-rs/commit/c8922fd4354e335bc138458916366d0d5e1c349a))
    - Making some options global so they can go anywhere on the cli ([`1ac197d`](https://github.com/fiberseq/fibertools-rs/commit/1ac197dde3e03289f2669829bfb2e1c1d651774b))
    - Making some options global so they can go anywhere on the cli ([`1f68cbe`](https://github.com/fiberseq/fibertools-rs/commit/1f68cbe5544c8eb0ac63f8aeada140bf56364dac))
    - Allow bam reading and writing from stdin for m6a-predict ([`c977510`](https://github.com/fiberseq/fibertools-rs/commit/c9775102d429f0e6ce0afae9219107b3f78aa4f0))
    - Add warning for missing hifi kinetics ([`21fe212`](https://github.com/fiberseq/fibertools-rs/commit/21fe212d8c51faffc74f7806cfc2f710f142aeb4))
    - Update install ([`3890f37`](https://github.com/fiberseq/fibertools-rs/commit/3890f379a35303edf2f51ebbcafb7863916f819d))
    - Update install ([`e4d0aad`](https://github.com/fiberseq/fibertools-rs/commit/e4d0aad2f8aeb273c9660e7150a366d2384d03a4))
</details>

## v0.0.8 (2022-11-20)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 68 commits contributed to the release over the course of 6 calendar days.
 - 11 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Prep for release ([`5a79467`](https://github.com/fiberseq/fibertools-rs/commit/5a7946757ea13b261e65bb0a5298442674149a62))
    - Prep for release ([`95c34ca`](https://github.com/fiberseq/fibertools-rs/commit/95c34ca6cc7bbffc6ab1637593c31e7b68fb4308))
    - Prep for release ([`f387eb4`](https://github.com/fiberseq/fibertools-rs/commit/f387eb432f8d22bf020c018396524cff791627fe))
    - Prep for release ([`7057bdf`](https://github.com/fiberseq/fibertools-rs/commit/7057bdf64fc40a3b345cef80f39f331c29e7b69b))
    - Check for RG ([`2a32603`](https://github.com/fiberseq/fibertools-rs/commit/2a32603dc7caaa3928de4a4fc84c29341fe8a49d))
    - Typo ([`812dff2`](https://github.com/fiberseq/fibertools-rs/commit/812dff2de3d75378d6e940aa4cff651362387209))
    - War messages ([`7ff93b0`](https://github.com/fiberseq/fibertools-rs/commit/7ff93b0420fdb67e3bd814f431f365b6ee203a20))
    - Updating docs for predict-m6a ([`41ac119`](https://github.com/fiberseq/fibertools-rs/commit/41ac119e82923f5148390c6c96b0fd6cb0a5eac7))
    - Updating docs for predict-m6a ([`d9aecda`](https://github.com/fiberseq/fibertools-rs/commit/d9aecdabad3fa75bfe3623e674830222c3ccaece))
    - Updating docs for predict-m6a ([`18213df`](https://github.com/fiberseq/fibertools-rs/commit/18213dff5bcb78c70d2b8ed3df96535378cafdf6))
    - Add predict to docs ([`4ccb774`](https://github.com/fiberseq/fibertools-rs/commit/4ccb7743d4aef11c2818a66d73e400982d3338e8))
    - Secondary check ([`8e34095`](https://github.com/fiberseq/fibertools-rs/commit/8e34095bb3be8642be75fec0e656bec888fce659))
    - Cleaning up ([`1ee535d`](https://github.com/fiberseq/fibertools-rs/commit/1ee535dd6212baee988b5dc15d64e3f5d84e6e64))
    - Filter for alignment types ([`d504543`](https://github.com/fiberseq/fibertools-rs/commit/d5045431cd4380267ed50aa0ea3e01b1503fe1f5))
    - Clippy ([`3c4088d`](https://github.com/fiberseq/fibertools-rs/commit/3c4088d449aec818cb0d00c04e72170d1a7ab7e9))
    - Rev qual? ([`603adb4`](https://github.com/fiberseq/fibertools-rs/commit/603adb48cdbe74b83e53f61182b0ff6e269f2714))
    - Rename ([`9ad12a7`](https://github.com/fiberseq/fibertools-rs/commit/9ad12a7ee31499b1a7dd91e579fcf4d000d29003))
    - Remove filter ([`0b58b77`](https://github.com/fiberseq/fibertools-rs/commit/0b58b77fe7f56680b2fbbb8edd61a4b8f8ac155f))
    - Updating basemod class ([`b86d481`](https://github.com/fiberseq/fibertools-rs/commit/b86d481431764a3c533f692f582524bd59674bbe))
    - New bigger ML models ([`00257fa`](https://github.com/fiberseq/fibertools-rs/commit/00257fa814c7a3794701e877863e0040f4ac7094))
    - I can now write ML MM tags using the BaseMods data structure ([`e0d546a`](https://github.com/fiberseq/fibertools-rs/commit/e0d546acdaa74d1613dc730dbf1da2093e16bb1f))
    - Filter mapped reads out of the process since it doesnt work right now ([`657b72a`](https://github.com/fiberseq/fibertools-rs/commit/657b72afbb7b5901c13e84d03048ca3536539c3a))
    - Check for missing liftovers, fiz ([`38a7dc3`](https://github.com/fiberseq/fibertools-rs/commit/38a7dc34202d11a63890bd4fcb5b297a1a20a6ba))
    - Check for missing liftovers, fiz ([`e91f54a`](https://github.com/fiberseq/fibertools-rs/commit/e91f54ad68a74162144e606a5e0116b1ac05673e))
    - Check for missing liftovers ([`0ce2c08`](https://github.com/fiberseq/fibertools-rs/commit/0ce2c0816ef9884f1b60091bc4db18bd5c59a0bf))
    - Fix prediction for mapped data? ([`ec6ec5a`](https://github.com/fiberseq/fibertools-rs/commit/ec6ec5aa91610df6bdb6dfd80a9f16eaac8b63c9))
    - Fix for real RC values with qual ([`db2aad9`](https://github.com/fiberseq/fibertools-rs/commit/db2aad93f753c078e1be7a6c7233a038873e5c7f))
    - Clean unused files and rename things in cnn ([`a340879`](https://github.com/fiberseq/fibertools-rs/commit/a340879415fab6ffb4cca167aa53466840a3b671))
    - Clean unused files and rename things in cnn ([`ce6335f`](https://github.com/fiberseq/fibertools-rs/commit/ce6335ffd2a8fad8af11d84ee4e3ea579cf12d35))
    - Clean unused files and rename things ([`4e7c710`](https://github.com/fiberseq/fibertools-rs/commit/4e7c710ffa395e9773cd3411bf8d003758b314da))
    - Remove filter ([`812a93e`](https://github.com/fiberseq/fibertools-rs/commit/812a93e3718ed719cddd3febade167b8248c4138))
    - Remove filter ([`bc4b80c`](https://github.com/fiberseq/fibertools-rs/commit/bc4b80c77fa7ce0d212e64abe7e5f8596dc11ff8))
    - Remove filter ([`9a7f5ff`](https://github.com/fiberseq/fibertools-rs/commit/9a7f5ff577f8c87cc7c5fe22f94b1de31320142d))
    - Remove filter ([`fc175ee`](https://github.com/fiberseq/fibertools-rs/commit/fc175eee6f26ea4040a25d1582998bf105250abe))
    - Remove filter ([`30a3119`](https://github.com/fiberseq/fibertools-rs/commit/30a31192b3a329405d91f062d6cd90f64064c522))
    - Remove filter ([`910014c`](https://github.com/fiberseq/fibertools-rs/commit/910014c15d3b9b4f1db027547e7384e4d7ba901c))
    - Remove filter ([`fbaf936`](https://github.com/fiberseq/fibertools-rs/commit/fbaf936eea3b965c0c10ddd187508aee90578cd0))
    - Remove filter ([`c7ab03d`](https://github.com/fiberseq/fibertools-rs/commit/c7ab03df80d2e2c6779fca75d14734894ca34dc8))
    - Remove filter ([`ce61ed5`](https://github.com/fiberseq/fibertools-rs/commit/ce61ed51d07f912a9bb17ac986c7ce14c60787b3))
    - Remove filter ([`e5df2bf`](https://github.com/fiberseq/fibertools-rs/commit/e5df2bf8fce36c6b5b6a294bb90e6b0fa2fe64ac))
    - Remove filter ([`b52f7e7`](https://github.com/fiberseq/fibertools-rs/commit/b52f7e7097374c29200eb2aa38430931103a3b8d))
    - Remove filter ([`e51d1c8`](https://github.com/fiberseq/fibertools-rs/commit/e51d1c8c988bc0ad5562bea38e817345ef39c5bc))
    - Remove filter ([`df4ef9a`](https://github.com/fiberseq/fibertools-rs/commit/df4ef9a4e9c9fb49e822c0dd8d8039f29fa9c9b2))
    - Remove filter ([`d8dcf92`](https://github.com/fiberseq/fibertools-rs/commit/d8dcf929ad5d4d87b390740651fbc9cef2cf94ef))
    - Remove filter ([`88a02c8`](https://github.com/fiberseq/fibertools-rs/commit/88a02c8b6e416ac586bb544cfda9c0902992a0f4))
    - Remove filter ([`a28e051`](https://github.com/fiberseq/fibertools-rs/commit/a28e0515f52fc0d5145217700d075e73118dd4bd))
    - Remove filter ([`9b2ea36`](https://github.com/fiberseq/fibertools-rs/commit/9b2ea36b8d0ab60e57230c8abe47a6aac4eeda00))
    - Remove filter ([`44ca3c7`](https://github.com/fiberseq/fibertools-rs/commit/44ca3c7d6c7369a6ad46dd954775be45624c4e41))
    - Remove filter ([`308f881`](https://github.com/fiberseq/fibertools-rs/commit/308f88133df22b1e0d60df1eb232d3894e46e1c2))
    - Assert ([`333494c`](https://github.com/fiberseq/fibertools-rs/commit/333494c30b0315e89ff249fecd87c1ede5ae8672))
    - Making a combined output for pos, ref_pos, and qual ([`3bbc7a1`](https://github.com/fiberseq/fibertools-rs/commit/3bbc7a16d8ff66fc2d53a41d08badcf3c1ab460c))
    - Making a combined output for pos, ref_pos, and qual ([`39e5d3a`](https://github.com/fiberseq/fibertools-rs/commit/39e5d3aee3a00d5b6dd8565ef2fdc3134483ee19))
    - Making a combined output for pos, ref_pos, and qual ([`537e38e`](https://github.com/fiberseq/fibertools-rs/commit/537e38e4f5c6b9fdd6b81a4f49854cbc57e2f7c6))
    - Making a combined output for pos, ref_pos, and qual ([`1c4c5b7`](https://github.com/fiberseq/fibertools-rs/commit/1c4c5b72d7096b975089cce4c6bff9d5927997dc))
    - Making a combined output for pos, ref_pos, and qual ([`34d3a82`](https://github.com/fiberseq/fibertools-rs/commit/34d3a8285e9b7cc839c012b7a07647fc5c68da52))
    - Making a combined output for pos, ref_pos, and qual ([`8c09829`](https://github.com/fiberseq/fibertools-rs/commit/8c09829e4de4763f12fa0c9cf1409227c9b5514d))
    - Another try ([`1d4d939`](https://github.com/fiberseq/fibertools-rs/commit/1d4d93971b17be4c11b6ab43ce061c479f4ecfe2))
    - Clear earlier ([`7d262ef`](https://github.com/fiberseq/fibertools-rs/commit/7d262efb1223cc7110c34619b973ace0c20d82d3))
    - Here ([`9d1a815`](https://github.com/fiberseq/fibertools-rs/commit/9d1a8155b189a5ef66331a2f6e39c0a4e1904a83))
    - Removing remove ([`488d5d5`](https://github.com/fiberseq/fibertools-rs/commit/488d5d5303d363d69f99de060ac010414242c343))
    - Moving basemods ([`ac2ecca`](https://github.com/fiberseq/fibertools-rs/commit/ac2eccae144d150a192ee85b41958c300e289412))
    - Fix qual? ([`00828e1`](https://github.com/fiberseq/fibertools-rs/commit/00828e11f51639608b97cabe5c00b6f54f608a61))
    - Fix qual? ([`8bd2e41`](https://github.com/fiberseq/fibertools-rs/commit/8bd2e41723aa2005495164ca0431193acf80d7e0))
    - Fix qual? ([`d69e61d`](https://github.com/fiberseq/fibertools-rs/commit/d69e61d0544d08b6d57925646f99ef26b73f627d))
    - Add missing 255 ([`048eb6d`](https://github.com/fiberseq/fibertools-rs/commit/048eb6dbc9a686cfb7030f9425963d275e5c655d))
    - Adding ability to extract full probablities from mp tag if it is there and requested ([`e981381`](https://github.com/fiberseq/fibertools-rs/commit/e981381b76cb8fb2a1265fef3cce3ae18353c2ef))
    - Optionally add full float ml prediction to bam ([`3b8f463`](https://github.com/fiberseq/fibertools-rs/commit/3b8f463f09ad8878c1040f7909983a959c341f63))
    - Add model for 2.2 chemsitry with xgb ([`8ae2a96`](https://github.com/fiberseq/fibertools-rs/commit/8ae2a96c5b63965c91eb8ebad679184effbc23f4))
</details>

## v0.0.7 (2022-11-08)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 105 commits contributed to the release over the course of 52 calendar days.
 - 54 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Check for unmapped ([`300ca87`](https://github.com/fiberseq/fibertools-rs/commit/300ca8748d63aeddbbd2bf68e13488b0dee5cb05))
    - Refactor xgb into its own file ([`1f38ace`](https://github.com/fiberseq/fibertools-rs/commit/1f38aced0b2c59d0e9135474f0d737297d8c092d))
    - Move cnn into a feature and its own file, and compile xgb predict by deafult ([`e4c3be2`](https://github.com/fiberseq/fibertools-rs/commit/e4c3be2c36feb5ee5c2c821a8ac085b2352e9be0))
    - Move cnn into a feature and its own file, and compile xgb predict by deafult ([`005327b`](https://github.com/fiberseq/fibertools-rs/commit/005327b5c878d850b813c4297b432097ed66a5ad))
    - Hide cnn predictions behind a feature flag, aka not compiled by default ([`8a811ca`](https://github.com/fiberseq/fibertools-rs/commit/8a811ca7780c8721847746157d046c681bb06f4a))
    - Updating docs with new help messages ([`a0dfed9`](https://github.com/fiberseq/fibertools-rs/commit/a0dfed979ecbdcd131a854fb6860224e6f81805a))
    - Adding seq remove option and quality addition option for extract-alll ([`6490336`](https://github.com/fiberseq/fibertools-rs/commit/6490336b9cd1b88f70725239098b35db8580e122))
    - Log scale scores ([`45ed225`](https://github.com/fiberseq/fibertools-rs/commit/45ed225ea3e67f2adf34f7bd3e470558c3f7f3cd))
    - Log scale scores ([`e1550ea`](https://github.com/fiberseq/fibertools-rs/commit/e1550eae9d4bef31f951783574fbe38b4d87c17e))
    - Log scale scores ([`e6775b9`](https://github.com/fiberseq/fibertools-rs/commit/e6775b94933c69e84cf205185beed39384741228))
    - Log scale scores ([`99c2307`](https://github.com/fiberseq/fibertools-rs/commit/99c23070b95f706c150d4eca8e45ee32cffc1e1d))
    - Log scale scores ([`a219299`](https://github.com/fiberseq/fibertools-rs/commit/a219299b53a08cf9782d6000622b77966b67479f))
    - Log scale scores ([`7fe2f64`](https://github.com/fiberseq/fibertools-rs/commit/7fe2f64ca0a3f9dfa1e0d6735081edeec719d477))
    - Log scale scores ([`fa3fe52`](https://github.com/fiberseq/fibertools-rs/commit/fa3fe52ebca70a948707c6ad6a510399d08edcef))
    - Log scale scores ([`d0c765f`](https://github.com/fiberseq/fibertools-rs/commit/d0c765f64ea07fcc6d5c7949596fcf1f26fdb699))
    - Log scale scores ([`7599a95`](https://github.com/fiberseq/fibertools-rs/commit/7599a957f9a2a55b784eca8d1d3fd1e7e6b424eb))
    - Log scale scores ([`cd19887`](https://github.com/fiberseq/fibertools-rs/commit/cd1988728a79595fdee65efa8eeae0accc5195d4))
    - Log scale scores ([`55c9d76`](https://github.com/fiberseq/fibertools-rs/commit/55c9d76d75630b5e72669cdd0a51177ab482c39c))
    - Log scale scores ([`7302915`](https://github.com/fiberseq/fibertools-rs/commit/73029154ad45b76b64833bc52ef2e4035bd68120))
    - Log scale scores ([`671fd39`](https://github.com/fiberseq/fibertools-rs/commit/671fd39bf67b32071abec8fd8610cd7da48de4c7))
    - Log scale scores ([`bace110`](https://github.com/fiberseq/fibertools-rs/commit/bace110f61504ce7e50e7e9507a045b7c5303638))
    - Log scale scores ([`b16e97a`](https://github.com/fiberseq/fibertools-rs/commit/b16e97a859ba0dcd724c47b4f99e941416e411b7))
    - Log scale scores ([`b896405`](https://github.com/fiberseq/fibertools-rs/commit/b896405d638955cfc173c946fe1792b76d82918a))
    - Log scale scores ([`fa6080e`](https://github.com/fiberseq/fibertools-rs/commit/fa6080ec45281967d1e41f88ff6dd1f516f45937))
    - Log scale scores ([`c5a6dc9`](https://github.com/fiberseq/fibertools-rs/commit/c5a6dc9360a8b7edffd68aa028cd80989e9dcf5d))
    - Log scale scores ([`85cdc97`](https://github.com/fiberseq/fibertools-rs/commit/85cdc97adb7c1901d06e44baff256bdee9106c88))
    - Log scale scores ([`682428c`](https://github.com/fiberseq/fibertools-rs/commit/682428cdb60717827b4364df2590cdfb62f59bf4))
    - Log scale scores ([`304cd12`](https://github.com/fiberseq/fibertools-rs/commit/304cd124dc39db3b7a543ba4cc3c11c92026ffb8))
    - Log scale scores ([`b1b07b6`](https://github.com/fiberseq/fibertools-rs/commit/b1b07b6e6d0fc7467200a711b7af9e3ccb339b6e))
    - Log scale scores ([`0f8605e`](https://github.com/fiberseq/fibertools-rs/commit/0f8605e86c329ac9e7d6232dee9365eb872d69e8))
    - Log scale scores ([`2d57517`](https://github.com/fiberseq/fibertools-rs/commit/2d57517b48f0a098a628bee0b065fabb1d8561ce))
    - Log scale scores ([`520405a`](https://github.com/fiberseq/fibertools-rs/commit/520405a4488297abe9223252ccd40419dd1f711d))
    - Log scale scores ([`3bdd3dc`](https://github.com/fiberseq/fibertools-rs/commit/3bdd3dc6ba3f3ea6276e24919f8eaf2791f6410f))
    - Log scale scores ([`0b1e12a`](https://github.com/fiberseq/fibertools-rs/commit/0b1e12a57d2f7cbd7882d750b9984e8dd9153612))
    - Log scale scores ([`b1cec34`](https://github.com/fiberseq/fibertools-rs/commit/b1cec34ec17a6ac686c12c52406fb4d6ad11a3ad))
    - Log scale scores ([`84df281`](https://github.com/fiberseq/fibertools-rs/commit/84df281b53031c3f6e4b54a5b3dbfb22857f50d3))
    - Log scale scores ([`9f9271e`](https://github.com/fiberseq/fibertools-rs/commit/9f9271e92f916b64c7cb69f9069a81d3c56cac87))
    - Log scale scores ([`8bb047c`](https://github.com/fiberseq/fibertools-rs/commit/8bb047ce971b1ad34d58d501269e8dae3e66c113))
    - Log scale scores ([`707cada`](https://github.com/fiberseq/fibertools-rs/commit/707cada84d493e8b637a155d871f109dd93631d3))
    - Log scale scores ([`368cf43`](https://github.com/fiberseq/fibertools-rs/commit/368cf438dd3e3cc59f2d5602fead89ef990cba83))
    - Log scale scores ([`efed35f`](https://github.com/fiberseq/fibertools-rs/commit/efed35f1e1e02701bec4becc2130671de27543ee))
    - Log scale scores ([`897cc0e`](https://github.com/fiberseq/fibertools-rs/commit/897cc0e9c917a7d4b6144322cea065b567983f91))
    - Log scale scores ([`0534a44`](https://github.com/fiberseq/fibertools-rs/commit/0534a449e3b1cfcbc35c0649f08defaef4b2b44a))
    - Log scale scores ([`2cc55d6`](https://github.com/fiberseq/fibertools-rs/commit/2cc55d6e1089fdf689bb43d886043305772b1d6b))
    - Log scale scores ([`9bc80b0`](https://github.com/fiberseq/fibertools-rs/commit/9bc80b0231095c1f0d5391c13724202ffd5a840d))
    - Log scale scores ([`e8ef551`](https://github.com/fiberseq/fibertools-rs/commit/e8ef55148364e0042ef70da970a9007553e0ff1d))
    - Log scale scores ([`ee5645d`](https://github.com/fiberseq/fibertools-rs/commit/ee5645dd2e35ffd86958bf59a0914367457a4a01))
    - Log scale scores ([`44a4a5d`](https://github.com/fiberseq/fibertools-rs/commit/44a4a5d17fcb1e410c70f1bbabe038f74c1474d3))
    - Log scale scores ([`5a02c02`](https://github.com/fiberseq/fibertools-rs/commit/5a02c0287b28ef3467514f94542e606b997e902e))
    - Log scale scores ([`7ef8d25`](https://github.com/fiberseq/fibertools-rs/commit/7ef8d256e62f9b541630fc87acc757ec5f4fd86f))
    - Log scale scores ([`bfb113f`](https://github.com/fiberseq/fibertools-rs/commit/bfb113fcf6b4fabaddb2a5eccca68a400f49dacd))
    - Log scale scores ([`686e6c4`](https://github.com/fiberseq/fibertools-rs/commit/686e6c4cb1af841c83caa3249148e644d3133dfd))
    - Log scale scores ([`6d9c943`](https://github.com/fiberseq/fibertools-rs/commit/6d9c94345b347b1b77e8b29505552678be9fe088))
    - Log scale scores ([`b0d8f3c`](https://github.com/fiberseq/fibertools-rs/commit/b0d8f3cac7ef61e2e4ee6092ea7103ae1e4f946c))
    - Log scale scores ([`79078ab`](https://github.com/fiberseq/fibertools-rs/commit/79078abbf5a019c987dcff0b774c3ef89047cb34))
    - Log scale scores ([`9049c72`](https://github.com/fiberseq/fibertools-rs/commit/9049c724f8f2d55e800b03cc078986672a794078))
    - Log scale scores ([`add7ba8`](https://github.com/fiberseq/fibertools-rs/commit/add7ba8711c0d77abcba1b2fad2edfa32bf28f02))
    - Log scale scores ([`fe90219`](https://github.com/fiberseq/fibertools-rs/commit/fe90219412724b1f6f8181f58a15bd136ccf9661))
    - Log scale scores ([`83cae68`](https://github.com/fiberseq/fibertools-rs/commit/83cae683ca71556c2fe6a71544a196451d56077e))
    - Log scale scores ([`523253c`](https://github.com/fiberseq/fibertools-rs/commit/523253c30aeff332d69cc20abcaa96acff8dbe88))
    - Log scale scores ([`490d1d2`](https://github.com/fiberseq/fibertools-rs/commit/490d1d2a00ed9c7cf1e5eeeb68cf27ee3822fd24))
    - Log scale scores ([`f5015d3`](https://github.com/fiberseq/fibertools-rs/commit/f5015d3490a20bd0aa0d7fe89f34ac4933a133d9))
    - Log scale scores ([`02c90f8`](https://github.com/fiberseq/fibertools-rs/commit/02c90f8007c33ffe6ba1995113c84ff5954a9fa6))
    - Log scale scores ([`aa39635`](https://github.com/fiberseq/fibertools-rs/commit/aa396350b08fc77e9a63d7997eecbe73ac063ddc))
    - Log scale scores ([`595907a`](https://github.com/fiberseq/fibertools-rs/commit/595907a0ccf59e34cc35ccdd0a824afacd948a72))
    - Log scale scores ([`25d2d20`](https://github.com/fiberseq/fibertools-rs/commit/25d2d205ce22f2376332c68124c0fcce5245ea23))
    - Log scale scores ([`1325001`](https://github.com/fiberseq/fibertools-rs/commit/132500194e3e2984e1496a09be87602bd4366a22))
    - Log scale scores ([`2319670`](https://github.com/fiberseq/fibertools-rs/commit/2319670bed73400af75c058248616a8719069c7b))
    - Log scale scores ([`6d366ce`](https://github.com/fiberseq/fibertools-rs/commit/6d366cea42b8c9be9c0c7e2c7f21714d8f381aa0))
    - Log scale scores ([`c1c5d7b`](https://github.com/fiberseq/fibertools-rs/commit/c1c5d7b64bfec0ab873c3969a107ff9ec11446db))
    - Fix hotone encodeing? ([`13f3e19`](https://github.com/fiberseq/fibertools-rs/commit/13f3e19f462a291b94496f154ee2a243421226e7))
    - Fix hotone encodeing? ([`549d7f1`](https://github.com/fiberseq/fibertools-rs/commit/549d7f1acb2914ad204e82e1953b57d9604450b6))
    - Fix offest in leading count? ([`5d6ffe2`](https://github.com/fiberseq/fibertools-rs/commit/5d6ffe2390336d061e2490102ffff7cf8b095db0))
    - Fix PW and IPW in the wrong direction on - strand ([`9d68bb2`](https://github.com/fiberseq/fibertools-rs/commit/9d68bb2ffbc5f4b449b50ec68e974c46f8c2afaa))
    - Update ([`f3e742c`](https://github.com/fiberseq/fibertools-rs/commit/f3e742c48283f39c7216e3b4ef8c5158e0094408))
    - Boost ([`ca24cfb`](https://github.com/fiberseq/fibertools-rs/commit/ca24cfb341d24cf1f5161a9a543a05f44654e687))
    - Cnn model update ([`9466260`](https://github.com/fiberseq/fibertools-rs/commit/946626075d8782e9758c1c1646243d42b2fc4db4))
    - Cnn model update ([`e235077`](https://github.com/fiberseq/fibertools-rs/commit/e235077bfb572f790b39a8b7a8a98f1d59da87e9))
    - Cnn model update ([`faeb253`](https://github.com/fiberseq/fibertools-rs/commit/faeb2536fb4defaf2cef55b20daea1e80162a171))
    - Cnn model update ([`b434fb1`](https://github.com/fiberseq/fibertools-rs/commit/b434fb16a6a91a237eb203a36f865db356a312ae))
    - Cnn model update ([`ba23e9b`](https://github.com/fiberseq/fibertools-rs/commit/ba23e9bef2473f04af3a08b76bac60030a3189b7))
    - Revetn cnn model ([`c6d907e`](https://github.com/fiberseq/fibertools-rs/commit/c6d907ece0a050663116912fb0cd4c74f39b7fe4))
    - New cnn model ([`cc7cdd5`](https://github.com/fiberseq/fibertools-rs/commit/cc7cdd518bf0acac682e736de6e7354e028e47bc))
    - Adding m6a qual to all output ([`e6bbdc8`](https://github.com/fiberseq/fibertools-rs/commit/e6bbdc85912ac2f471ea745cb4797cc1e8e204f8))
    - Adding a new m6a call that also gets quals ([`93a005f`](https://github.com/fiberseq/fibertools-rs/commit/93a005fa1813b9571010429472bf7afa6bf6112b))
    - New test ([`b74871c`](https://github.com/fiberseq/fibertools-rs/commit/b74871c1a344a4b1d4c0e77fca98bcc8adfccc92))
    - New test ([`17dbc9d`](https://github.com/fiberseq/fibertools-rs/commit/17dbc9dcd430d0bd90635c240309f31a05d0595b))
    - Larged model ([`f10c4d3`](https://github.com/fiberseq/fibertools-rs/commit/f10c4d3039996a5cce096ca62dda3f693ea3d7db))
    - Adding comple flags for tch? ([`b81248a`](https://github.com/fiberseq/fibertools-rs/commit/b81248aec0e5fcf7235f34c2001577b5de2b6456))
    - Adding comple flags for tch? ([`a6c324f`](https://github.com/fiberseq/fibertools-rs/commit/a6c324ff3bbcb9533a70b24ddafb7b38607054f4))
    - Adding a cnn flag ([`1184418`](https://github.com/fiberseq/fibertools-rs/commit/1184418abfe60ba9a544843dee371e4cbfae47f2))
    - Adding the cnn model? ([`f085bac`](https://github.com/fiberseq/fibertools-rs/commit/f085bac197e84e3a38deb8c8902d0b6ee039d64a))
    - Adding out threads ([`da1f980`](https://github.com/fiberseq/fibertools-rs/commit/da1f98012e79b7b07a29f5cae0ccaf30541645b0))
    - Remove tmep ([`458a062`](https://github.com/fiberseq/fibertools-rs/commit/458a06204e6b6e4f8bd4811b5538a8b8cfa6fe77))
    - Log ([`6e8f9e6`](https://github.com/fiberseq/fibertools-rs/commit/6e8f9e6b73c831fed0efe530029bd98a75cab47d))
    - Clippy ([`5670402`](https://github.com/fiberseq/fibertools-rs/commit/56704026d83a00fdae7eaad4213ab850c65e336b))
    - Moved to a rust xgboost ([`0af3137`](https://github.com/fiberseq/fibertools-rs/commit/0af3137ee15d51afc025f7f422515ba2442f6d58))
    - Added a basic m6a caller that adds to the bam tags ([`71e5812`](https://github.com/fiberseq/fibertools-rs/commit/71e58128417c909190726b624e5e1d32eda665a6))
    - Clean ([`a723d61`](https://github.com/fiberseq/fibertools-rs/commit/a723d613358b6131185bb17a93b5fb4924e1dd4f))
    - Added a proof of principle xgboost application to predicting m6A ([`e70e7b9`](https://github.com/fiberseq/fibertools-rs/commit/e70e7b956c36b5e2323370b0e78de57566d39802))
    - Smaller ([`ab1ea68`](https://github.com/fiberseq/fibertools-rs/commit/ab1ea68bc5467c7a3a018a99b9a277bf755a9723))
    - Versions ([`121e7c9`](https://github.com/fiberseq/fibertools-rs/commit/121e7c9837b76e2f7b7c72937bfa679881f06f14))
    - Versions ([`72f7cdd`](https://github.com/fiberseq/fibertools-rs/commit/72f7cdded1903587ff2d839a4c4f91a1fd10290b))
    - Versions ([`184fe65`](https://github.com/fiberseq/fibertools-rs/commit/184fe650e8a3d7ea03b97fe2c9b11c58ef6f65a5))
    - Ignore ([`1a7a741`](https://github.com/fiberseq/fibertools-rs/commit/1a7a741f624dedfe4b3a5b32212f82db4e08ce9c))
</details>

## v0.0.6 (2022-09-15)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 13 commits contributed to the release over the course of 14 calendar days.
 - 14 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Makeing the 00.0.6 release ([`8847fe5`](https://github.com/fiberseq/fibertools-rs/commit/8847fe5e133b62eed66915ef5dfbcfb1dda2712b))
    - Updating readme and docs ([`bd22fb1`](https://github.com/fiberseq/fibertools-rs/commit/bd22fb1ccf68614e3b211eafa2f07081ce513985))
    - Updating readme and docs ([`b086baa`](https://github.com/fiberseq/fibertools-rs/commit/b086baa4644098c153672fd2fe5f255f80f8f902))
    - Adding rq tag to all ([`9ae61bb`](https://github.com/fiberseq/fibertools-rs/commit/9ae61bb338ff89da9c3e0aa556c91ec6a9cd95c3))
    - Todos ([`28775b6`](https://github.com/fiberseq/fibertools-rs/commit/28775b6bb9762847ec89ff5c85ce71b74ba24dcc))
    - Fixing CI ([`279dcc2`](https://github.com/fiberseq/fibertools-rs/commit/279dcc2f252e861cb25113ffdd842fdf4616d226))
    - Changing test file and expanding CI ([`44e4cb7`](https://github.com/fiberseq/fibertools-rs/commit/44e4cb7d56a12b4c4c1deb0b974db928481e3bb1))
    - Changing test file and expanding CI ([`cfe8e4c`](https://github.com/fiberseq/fibertools-rs/commit/cfe8e4c7145511f32e5f28fc71fbc6dfcd03152a))
    - Fixing CI ([`cc5efd2`](https://github.com/fiberseq/fibertools-rs/commit/cc5efd2b02ad3c2c5c11179d4c3ffaeb5554e671))
    - Adding CI badge ([`4bbeb45`](https://github.com/fiberseq/fibertools-rs/commit/4bbeb45e9c7ab54c2b761932891c961536b179af))
    - Adding CI ([`71b3ade`](https://github.com/fiberseq/fibertools-rs/commit/71b3adef427d601d1e96632b0403654490a2c029))
    - Adding CI ([`03efba5`](https://github.com/fiberseq/fibertools-rs/commit/03efba5e4b994f58f8dd8bfd3fe29e26a44b8925))
    - Adding CI ([`19b024d`](https://github.com/fiberseq/fibertools-rs/commit/19b024d887163ff66fe55342b737b4e5425ac7af))
</details>

## v0.0.6-alpha.1 (2022-08-31)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 11 commits contributed to the release over the course of 5 calendar days.
 - 7 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Making alpha version for cargo ([`0a68568`](https://github.com/fiberseq/fibertools-rs/commit/0a68568e0f097080e870c50239ca0779ef96d431))
    - Making alpha version for cargo ([`818c6e0`](https://github.com/fiberseq/fibertools-rs/commit/818c6e05694e85c1b5e0bc7bef4f65dc773cfd64))
    - Changeing writer to handle broken pipes, line with head, and removing some completed TODOs ([`074bcab`](https://github.com/fiberseq/fibertools-rs/commit/074bcab25a0874e81975042cb918d6cd8c9b5d60))
    - Changeing writer to handle broken pipes, line with head, and removing some completed TODOs ([`8c7e31b`](https://github.com/fiberseq/fibertools-rs/commit/8c7e31b5278d92136b12f7311a771c4ad3602179))
    - TODO ([`fa2c31c`](https://github.com/fiberseq/fibertools-rs/commit/fa2c31ce7f484d9551eac2a291c776331b148350))
    - Fixing a suspected off by one issue on the reverse strand ([`63840dc`](https://github.com/fiberseq/fibertools-rs/commit/63840dcf47913ff01b43ed2f21b516fcb7f863a0))
    - Adding total m6a, AT, and 5mc counts ([`32397e8`](https://github.com/fiberseq/fibertools-rs/commit/32397e8bab0d17e624c4c6a4c3454d342df04048))
    - Replace missing contig with a . ([`1ea4769`](https://github.com/fiberseq/fibertools-rs/commit/1ea476947c79dab1570fec1a13db3cd8b969d8ab))
    - Fixed segfault for unaligned reads trying to get reference info ([`ea05a59`](https://github.com/fiberseq/fibertools-rs/commit/ea05a59eb7c1ff485859bcd1e1ebccd554a46241))
    - Improve install instructions ([`cc302d4`](https://github.com/fiberseq/fibertools-rs/commit/cc302d4f8baa88168d97fc56a3822998785c4b5c))
    - Checking for unmapped in ranged liftover ([`52eae4b`](https://github.com/fiberseq/fibertools-rs/commit/52eae4bd1785756b1e23e723c7691714c0f4a5bf))
</details>

## v0.0.5 (2022-08-23)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 29 commits contributed to the release over the course of 23 calendar days.
 - 23 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Bump ([`b2c5663`](https://github.com/fiberseq/fibertools-rs/commit/b2c5663bda4b879fc635386ad858817d80a2d961))
    - Adding links to automated help pages to help keep uptodate. Closes #1 ([`2ee0542`](https://github.com/fiberseq/fibertools-rs/commit/2ee0542ba6fdd68b0a12014808554ed96c0ae268))
    - Made a large # of changes so that outputs are consistently in the orientation of the read stored in the bam rather than the original orientation ([`7c6f773`](https://github.com/fiberseq/fibertools-rs/commit/7c6f773d50a5d60a7dea4de7895ee7381115a90b))
    - Fixed bug for missing reverse strand outputs ([`318aac5`](https://github.com/fiberseq/fibertools-rs/commit/318aac56263db53e570b91db0ef6bc77e80d8357))
    - Removing filter to keep the -1 ([`38c8a77`](https://github.com/fiberseq/fibertools-rs/commit/38c8a77f6e53c218c923bfafe2e3327839411add))
    - Renaming output to 5mC ([`23d4526`](https://github.com/fiberseq/fibertools-rs/commit/23d45268bc3df175e76e86b46c7f1f7fc3735bad))
    - Adding a -1 value for positions that cannot be lifted over in extract all ([`b50df1b`](https://github.com/fiberseq/fibertools-rs/commit/b50df1b064ba42ac591cee72182dc0d910f52fad))
    - Adding . for missing feilds ([`6f6ee28`](https://github.com/fiberseq/fibertools-rs/commit/6f6ee28a163fcd9e70f149f92b590e2fa64366b5))
    - Remove tab at end of line ([`d52b388`](https://github.com/fiberseq/fibertools-rs/commit/d52b388570f5bc9e3e59af4b23ffe972afb7dce5))
    - Removing duplicate headers ([`81de168`](https://github.com/fiberseq/fibertools-rs/commit/81de168625747b52c31e692c519a591b8c3ca964))
    - Adding img ([`e219b79`](https://github.com/fiberseq/fibertools-rs/commit/e219b7909f8baebaad266872135cf77867b60d1e))
    - Adding comment char for all header ([`07be127`](https://github.com/fiberseq/fibertools-rs/commit/07be127be30e3a229bdd0a5d83baeab947618e7e))
    - Adding column names for --all ([`1f00737`](https://github.com/fiberseq/fibertools-rs/commit/1f007377ef92cb0c7cd91c787d61442cf5b02a88))
    - Try again ([`f08170f`](https://github.com/fiberseq/fibertools-rs/commit/f08170f82b0a57345318bc6b3f1a3ba4e3dab7a6))
    - Try again ([`e336377`](https://github.com/fiberseq/fibertools-rs/commit/e3363778e1e2085a5c52660ce06db5b8bc95cec3))
    - Try again ([`20b2050`](https://github.com/fiberseq/fibertools-rs/commit/20b20501c15c90a71d8892e10f17809273c005da))
    - Try again ([`66a506f`](https://github.com/fiberseq/fibertools-rs/commit/66a506f6eb8da407a50a4e6589ffdf5a2bdb6a25))
    - Try again ([`7e7cad3`](https://github.com/fiberseq/fibertools-rs/commit/7e7cad3d8313fbd673226dab78c8a36163752332))
    - Try again ([`fe0a6b4`](https://github.com/fiberseq/fibertools-rs/commit/fe0a6b4d65ce14a003f34461f65acd809bee6b3a))
    - Try again ([`96f231a`](https://github.com/fiberseq/fibertools-rs/commit/96f231ac977f1f53c3adf12e4e4d8dd9d5ccc526))
    - Try again ([`01e9ede`](https://github.com/fiberseq/fibertools-rs/commit/01e9edebe6d9967a9dc7d70b4880f2ba77229b9e))
    - Try again ([`2f3a2fb`](https://github.com/fiberseq/fibertools-rs/commit/2f3a2fb433623e6fe7a18b79411ea9fe768f5b37))
    - Try again ([`c27f71d`](https://github.com/fiberseq/fibertools-rs/commit/c27f71d4a5be58a2416849b97f3b28c46d0bdcd1))
    - Try again ([`ae31d3c`](https://github.com/fiberseq/fibertools-rs/commit/ae31d3cd95a2fef2f160fc2e7cfc4df7a3d1b21a))
    - Try again ([`c65d806`](https://github.com/fiberseq/fibertools-rs/commit/c65d806784814b0b994fddff711d5300d154c927))
    - Try again ([`034e34a`](https://github.com/fiberseq/fibertools-rs/commit/034e34a387c0d534c65d2be376acf6ec4f243cc4))
    - Try again ([`17d6705`](https://github.com/fiberseq/fibertools-rs/commit/17d6705d5b673819fcf8fda60e737bc10be5bbdd))
    - Try again ([`fa6ee99`](https://github.com/fiberseq/fibertools-rs/commit/fa6ee99a5a00fcdd4df12abe0f06283451fa99f0))
    - Try again ([`26132de`](https://github.com/fiberseq/fibertools-rs/commit/26132de803b678c6981408e722a5292d05edc3b2))
</details>

## v0.0.4 (2022-07-31)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 31 commits contributed to the release over the course of 2 calendar days.
 - 2 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Try again ([`aba1393`](https://github.com/fiberseq/fibertools-rs/commit/aba1393448fa008ba61400fd2fce05560b06df4f))
    - Try again ([`b0d39db`](https://github.com/fiberseq/fibertools-rs/commit/b0d39dbdac9bda011983f2e185fab17f025f5c74))
    - Try again ([`d460c52`](https://github.com/fiberseq/fibertools-rs/commit/d460c52190a47fd06c90960a88e932cf074fb91f))
    - Try again ([`c6ea96b`](https://github.com/fiberseq/fibertools-rs/commit/c6ea96b20ac4906b827e626ed196bf5104d0c978))
    - Try again ([`ff04b4d`](https://github.com/fiberseq/fibertools-rs/commit/ff04b4de7fafcc07bf378fddb3af38ca35aa3fb0))
    - Try again ([`2ac517b`](https://github.com/fiberseq/fibertools-rs/commit/2ac517b8ecf997fc8b54b4946d546eabaf6798d6))
    - Try again ([`e2faa05`](https://github.com/fiberseq/fibertools-rs/commit/e2faa0543a49ef7b7c2abb6f69d7b5ecba5b0eb8))
    - Try again ([`f998631`](https://github.com/fiberseq/fibertools-rs/commit/f998631f15f30993e2f430dba6966fd247fa84fa))
    - Try again ([`a876050`](https://github.com/fiberseq/fibertools-rs/commit/a8760507fc547e026d6f47f977bb6ba032bb0e5b))
    - Try again ([`887fec5`](https://github.com/fiberseq/fibertools-rs/commit/887fec55bfeeeb6ff65cb8b83a4743b2b1e3eb13))
    - Try again ([`5b10de4`](https://github.com/fiberseq/fibertools-rs/commit/5b10de4241691f8940721c9040aeaa6c065e96ea))
    - Try again ([`52d7bd7`](https://github.com/fiberseq/fibertools-rs/commit/52d7bd714c8b0de7f3db5895d8ea00a7a54aaadf))
    - Try again ([`4e05b62`](https://github.com/fiberseq/fibertools-rs/commit/4e05b62480ff1e70df6d76639962a4085c44b693))
    - Try again ([`20d3e10`](https://github.com/fiberseq/fibertools-rs/commit/20d3e103b93baa19bc08ffb78a5e0d7a326664ef))
    - Try again ([`e7f890c`](https://github.com/fiberseq/fibertools-rs/commit/e7f890cb77b8e7b841a3f304ad4e047c91be0d5a))
    - Try again ([`cbd77e4`](https://github.com/fiberseq/fibertools-rs/commit/cbd77e473f4bf364e45759d93359d396fb4569f5))
    - Try again ([`4c88f74`](https://github.com/fiberseq/fibertools-rs/commit/4c88f74e74a14fbafc632663e02276c86bb126f4))
    - Try again ([`86c7202`](https://github.com/fiberseq/fibertools-rs/commit/86c720283166eaa6d6687e8b7cffcf5090d79eec))
    - Try again ([`b090f10`](https://github.com/fiberseq/fibertools-rs/commit/b090f10710645fbd7e023a145c852590be9438c5))
    - Try again ([`4b2e4f6`](https://github.com/fiberseq/fibertools-rs/commit/4b2e4f620ea8fa439e946494919d02b79c1f9d23))
    - Try again ([`b47f660`](https://github.com/fiberseq/fibertools-rs/commit/b47f660e7418121a29410e187270828f0b7bcf5f))
    - Try again ([`42aba01`](https://github.com/fiberseq/fibertools-rs/commit/42aba0163ab75b56b7744bbfc59a57949391369c))
    - Try again ([`47b266e`](https://github.com/fiberseq/fibertools-rs/commit/47b266ebba89d06c4d611dc3ae30f16067604d9f))
    - Try again ([`b21de47`](https://github.com/fiberseq/fibertools-rs/commit/b21de47fe60df83b92e4e76289e6f1e67434a749))
    - Try again ([`d394151`](https://github.com/fiberseq/fibertools-rs/commit/d394151025a78c8ca888d5bf41e7df7995499a28))
    - Try again ([`0b27603`](https://github.com/fiberseq/fibertools-rs/commit/0b2760336c609df0f06687b9f4c7e93c641e0cee))
    - Try again ([`c60bcc5`](https://github.com/fiberseq/fibertools-rs/commit/c60bcc5db5bab3d23afda04fb84e194dd7236d2a))
    - Try again ([`4b0ffca`](https://github.com/fiberseq/fibertools-rs/commit/4b0ffcafa9e84a6c0dc18444404a05d70f498963))
    - Try again ([`852a562`](https://github.com/fiberseq/fibertools-rs/commit/852a56290e03e45fce6e5964c9f6a7156c9df25a))
    - Try again ([`40cfdb0`](https://github.com/fiberseq/fibertools-rs/commit/40cfdb0d6a4d35e182e74110da35c86feef21aab))
    - Try again ([`bc3c10f`](https://github.com/fiberseq/fibertools-rs/commit/bc3c10fb8dd47ee27f799b49e6d0c215a1d32693))
</details>

## v0.0.3-alpha (2022-07-28)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 18 commits contributed to the release over the course of 1 calendar day.
 - 1 day passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Try again ([`f8d57b0`](https://github.com/fiberseq/fibertools-rs/commit/f8d57b0b7c24db853d0554a6552b1e54c44ea8d7))
    - Try again ([`a7b022e`](https://github.com/fiberseq/fibertools-rs/commit/a7b022ed10d2cf5db382b5ad3b0b3c56e398fe1e))
    - Try again ([`8f5e64a`](https://github.com/fiberseq/fibertools-rs/commit/8f5e64aab138e5ba6661736e02f512d354f60cc7))
    - Try again ([`17fb908`](https://github.com/fiberseq/fibertools-rs/commit/17fb908a5d40e5d43b8363f3d540b492a78c83ba))
    - Try again ([`e961e4c`](https://github.com/fiberseq/fibertools-rs/commit/e961e4c1043e7b9d32fe7f3f4f7a9a0051eb6874))
    - Try again ([`c2baedf`](https://github.com/fiberseq/fibertools-rs/commit/c2baedf4e8096dd70c24d0b2a7f412941ac37eeb))
    - Try again ([`a1c4a3c`](https://github.com/fiberseq/fibertools-rs/commit/a1c4a3cfbb6f41ee1636bda46388b7a43268a7f1))
    - Try again ([`0beb42e`](https://github.com/fiberseq/fibertools-rs/commit/0beb42edf4758812ae39f8f26907cef7872d136b))
    - Try again ([`010dbfd`](https://github.com/fiberseq/fibertools-rs/commit/010dbfde2f5dc91f6f63f16cf32384cdd76d3c91))
    - Try again ([`b006a78`](https://github.com/fiberseq/fibertools-rs/commit/b006a78163709eb3682b6ff408ad4f3f6d217e5c))
    - Try again ([`8b328bf`](https://github.com/fiberseq/fibertools-rs/commit/8b328bf3cf5b803c3b6f4f12ed2899b28918b174))
    - Try again ([`5001a71`](https://github.com/fiberseq/fibertools-rs/commit/5001a719cc12ccd03c7f02af70062ed15473977b))
    - Try again ([`c798640`](https://github.com/fiberseq/fibertools-rs/commit/c79864032038434534890db1d9686ad8566e8f2e))
    - Try again ([`e8e1a38`](https://github.com/fiberseq/fibertools-rs/commit/e8e1a38ec3ce708a85fb3e27e7b8b7352f4f93f1))
    - Try again ([`13ad33f`](https://github.com/fiberseq/fibertools-rs/commit/13ad33f58703bcc4e760e979f4684876c65b9d0a))
    - Try again ([`de3267d`](https://github.com/fiberseq/fibertools-rs/commit/de3267d6d2d5396a98fa238ba5fe8327363f9b30))
    - Try again ([`967533c`](https://github.com/fiberseq/fibertools-rs/commit/967533ce1a4c495ea89229962d1d30bb7c95156b))
    - Try again ([`a9ebd24`](https://github.com/fiberseq/fibertools-rs/commit/a9ebd24948d303391c82d761246d49e9867426ce))
</details>

## v0.0.2 (2022-07-27)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 6 commits contributed to the release.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Try again ([`157a2aa`](https://github.com/fiberseq/fibertools-rs/commit/157a2aafd1ce4ee25301d571f5e32d68e04dbf2c))
    - Try again ([`6ac49c6`](https://github.com/fiberseq/fibertools-rs/commit/6ac49c6ae3dc43fc283faebdfc1b7d20973d99f6))
    - Build# ([`7caf415`](https://github.com/fiberseq/fibertools-rs/commit/7caf4156dd4a385c275e679b8f1a88d021243d28))
    - Adding ([`bfafb3b`](https://github.com/fiberseq/fibertools-rs/commit/bfafb3b812b8b4672c2ff278253e0032b8603067))
    - Adding ([`5f03308`](https://github.com/fiberseq/fibertools-rs/commit/5f03308482d48c34f144b197aa8b935f21a820d4))
    - Set theme jekyll-theme-minimal ([`3c3b11f`](https://github.com/fiberseq/fibertools-rs/commit/3c3b11f7087a90f095d0d610732e919fdac46746))
</details>

## v0.0.1 (2022-07-27)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 35 commits contributed to the release over the course of 2 calendar days.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Adding ([`493de98`](https://github.com/fiberseq/fibertools-rs/commit/493de98db23e567607d496422ece6ac89bfb497b))
    - Adding ([`f0d4310`](https://github.com/fiberseq/fibertools-rs/commit/f0d43109d2101a9e415d95e1a6825bf62c74ef7b))
    - Adding ([`116af9d`](https://github.com/fiberseq/fibertools-rs/commit/116af9df584a21ff0e0b9634f2a8e08d64ddb204))
    - Adding ([`5e2dc3b`](https://github.com/fiberseq/fibertools-rs/commit/5e2dc3b62dd71d4c3c63d5ea06d1ff7a17fe9b66))
    - Adding ([`5bb45b3`](https://github.com/fiberseq/fibertools-rs/commit/5bb45b397f817ce14310fdbd953ea03c99085124))
    - Adding ([`81f6dfd`](https://github.com/fiberseq/fibertools-rs/commit/81f6dfd92a81976492fe2d0ec5eddc063bb3111d))
    - Adding ([`17caed7`](https://github.com/fiberseq/fibertools-rs/commit/17caed707b099304802f3195d95beb5dbbf23e8d))
    - Adding ([`131fb08`](https://github.com/fiberseq/fibertools-rs/commit/131fb0826cf7567395846d80db553ece267ba163))
    - Adding ([`f1a51e9`](https://github.com/fiberseq/fibertools-rs/commit/f1a51e91255515bf0e28cae8f9dc65e77f477e3a))
    - Adding ([`66951c7`](https://github.com/fiberseq/fibertools-rs/commit/66951c7b1f48425b97ab540e12ffca074ea94795))
    - Adding ([`7c77792`](https://github.com/fiberseq/fibertools-rs/commit/7c777921ca3bd7dd29a71dd1b9d4aa765d238180))
    - Adding ([`7ff38f9`](https://github.com/fiberseq/fibertools-rs/commit/7ff38f95c2a81811581b84706bc67974ff739929))
    - Adding ([`e9234cd`](https://github.com/fiberseq/fibertools-rs/commit/e9234cddbc97ca2de160c5f98ea677cfa6b1781a))
    - Adding ([`316bcbb`](https://github.com/fiberseq/fibertools-rs/commit/316bcbb578e6751c68ef1ed3300afdbe1426d2e2))
    - Adding ([`1cbf8ab`](https://github.com/fiberseq/fibertools-rs/commit/1cbf8ab50a585d3bbbdfe0517683d5271849b2f4))
    - Adding ([`02c9983`](https://github.com/fiberseq/fibertools-rs/commit/02c9983a5f55f57f93f2c739cfb4774218e82e03))
    - Adding ([`cb24e21`](https://github.com/fiberseq/fibertools-rs/commit/cb24e21a5f72a4e166be88c2f4bbe25d1ec87297))
    - Adding ([`9bf5ebc`](https://github.com/fiberseq/fibertools-rs/commit/9bf5ebc2930d4d9ba9b586f7294ffe63de185ec4))
    - Adding ([`a7f3707`](https://github.com/fiberseq/fibertools-rs/commit/a7f3707d0bf8aa3c0e55079c05c9a194b59486a9))
    - Adding ([`5984937`](https://github.com/fiberseq/fibertools-rs/commit/5984937945abcb7218f0f15fdeec9ba5a1f46045))
    - Adding ([`e21e345`](https://github.com/fiberseq/fibertools-rs/commit/e21e34559b3b176f664a8fe511284353e6763126))
    - Adding ([`d5fb531`](https://github.com/fiberseq/fibertools-rs/commit/d5fb531a6fba4466a08ed343e817ed939fbe7def))
    - Adding ([`dad97a2`](https://github.com/fiberseq/fibertools-rs/commit/dad97a222bfb05003b441b9adcc92dc45af6848d))
    - Adding ([`1384ff1`](https://github.com/fiberseq/fibertools-rs/commit/1384ff1ae8b3d5918ebf937c9c2971de2edae1e3))
    - Adding ([`1040d32`](https://github.com/fiberseq/fibertools-rs/commit/1040d32c250fec34c048b051d2350baabd3ca45a))
    - Adding ([`773be4e`](https://github.com/fiberseq/fibertools-rs/commit/773be4e7e183f4c26b0c7367f6f3828ec2251b2b))
    - Adding ([`66b72f8`](https://github.com/fiberseq/fibertools-rs/commit/66b72f8ce9993974740c0e273f360db305ea0719))
    - Adding ([`d3cbad4`](https://github.com/fiberseq/fibertools-rs/commit/d3cbad4d369f0466e3c45ac17e2555edbfcd3dd9))
    - Adding ([`7045b27`](https://github.com/fiberseq/fibertools-rs/commit/7045b27483d453934956b66453d02769196d1c2b))
    - Adding ([`5bbadf7`](https://github.com/fiberseq/fibertools-rs/commit/5bbadf78d149f487ad7e4e6a3bb0ad21635e8682))
    - Adding ([`cc72a38`](https://github.com/fiberseq/fibertools-rs/commit/cc72a38bf8ed015355eb192398de65b5a9b89dfd))
    - Adding ([`5a81075`](https://github.com/fiberseq/fibertools-rs/commit/5a810759ed4a2bd6f908c83d3d6348d0282bf5c0))
    - Adding ([`832dd42`](https://github.com/fiberseq/fibertools-rs/commit/832dd427a8d560aba24c50529c98fb584f0fa080))
    - Adding ([`477bc7e`](https://github.com/fiberseq/fibertools-rs/commit/477bc7e7323d923b6ff55bc003d28877103ff5b4))
    - First commit ([`13746b2`](https://github.com/fiberseq/fibertools-rs/commit/13746b24d9b5a323c625909c61a1edcfcb1514aa))
</details>


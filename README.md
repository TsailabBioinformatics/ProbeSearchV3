# ProbeSearchV3
[//]: <TODO: figure out the best way to structure and tranfer `indices/` and `data/` directories>

This app consists of a form, which allows the user to input a primer/probe/gRNA sequence and select a set of genome databases. Upon form submission, the app returns a result of the input read aligned against the selected genome. The alignment is done via bowtie2 software. The app is written in Nodejs, Expressjs, and Vue.js. 

### Requirements
- Node.js: https://nodejs.dev/learn/how-to-install-nodejs


### Instructions
1. Clone this repository. \
`git clone https://github.com/TsailabBioinformatics/ProbeSearchV3`

2. Travel to the parent project directory and install the required node packages. \
`npm install` \
Then change into the frontend `vue/` directory, and install its packages. \
`cd vue/` \
`npm install`

3. Make the directories `data/` and `indices/` directories to the parent directory. Refer to [Editing Backend section.](#editting-probesearchv3)

4. Run the app. \
`node app.js` 

---

### Understanding ProbeSearchV3
Here is description of full stack of ProbeSearchV3. 

**Directory Structure**
- `data/`: Fasta files are located in `data/`
- `bowtie/`: Bowtie2 executables are located in `bowtie/`  
- `scripts/`: Helper scripts for Bowtie2
- `indices/`: Bowtie2 uses indexed versions of the fasta files. Indices are located in `indices/`. 
  - In order to create an index, use `./bowtie/bowtie2-build <path_to_fasta> <index_name>`. Then move these indices to `indices/`
  - Indices are integral to bowtie2 alignment. `app.js` makes a call to `./bowtie/bowtie2 -x indices/<db> -k 30 -c <read>`. Here, `<db>` represents the user-selected database, and `<read>` represents the input read. 
- `vue/`: Frontend Vue files

**Data Flow** \
The following details how data moves in the app. We'll see how the client-side form submission *request* leads to the eventual alignment result *response*. 

1. Initial Input Processing \
  User input is first handled by `SequenceForm.vue`. Within `SequenceForm.vue`, these data are tracked as `read: ''` and `db: []`. `read` is a string because it represents a primer/probe/gRNA sequence, and `db` is an array because it represents a set of genome databases. When the user fills out the form, the variables are defined respetively. Then, when the user submits the form (e.g., clicks `Search`), the variables are propagated, or *emitted* (in Vue lingo), up to its parent component `App.vue`. `App.vue`, then processes the user input and makes a put request to the Express api, `app.js` (`axios.put('/', payload)`). 
  
2. API Functions \
  The put function first calls a child process: `./bowtie/bowtie2 -x indices/ + String(req.body.db) -k 30 -c + String(req.body.read) + --end-to-end --no-hd`. Note, `req.body.db` and `req.body.read` is the Express way of accessing the inputted db and read respectively. This child process initiates the bowtie2 alignment function, which returns a SAM file as `stdout`. Then, the SAM file, or `stdout` is parsed and molded into an alignment illustration. The illustration of the alignmened is thne returned to `App.vue`. 
  
3. Returning Alignment Result \
  `app.js` sends back data in the form of an illustration to `App.vue`. In turn, `App.vue` passes this data to a child component called `AlignmentResult.vue`. Once the `AlignmentResult` component(s) receive the data, they show on client side. 

**Deployment** \
TBD, but as of now, the app is deployed at:  http://aspendb.uga.edu/probesearch/v3/

#### Editting ProbeSearchV3
- **Frontend**
  - The frontend can be editted from within the `vue/` directory. 
  - The three components as of now are the parent, `App.vue`, and its children, `SequenceForm.vue` & `AlignmentResult.vue`. 
  - After editting the vue files, run `npm run build` from within the `vue/` directory. This command builds the static HTML, css, and js for the frontend and places them in `vue/dist/` directory. The express app, `app.js`, references these static files automatically, so that's all: edit then build. 
- **API**
  - The main app can be editted at `app.js`. It is written in Express.
- **Backend**
  - The backend consists of two mains parts: the fasta files & the indexed fasta files. In order to edit the backend, or implement more genomes, do the following:
  1. Transfer whatever fasta files you wish to the parent directory.
  2. Build the indices: `./bowtie/bowtie2-build <path_to_fasta> <index_name>`. Name the index (second parameter) the same name as you display the database on the frontend.
  3. Move the raw fasta file to the `data/` directory, and move the indexed fasta files to the `indices/` directory.
  4. Go the Express `app.js`, and add to `db_dictionary`, which facilitates a dictionary between the frontend name of the database and the location of the fasta file.

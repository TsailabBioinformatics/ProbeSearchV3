<script>
  export default {
    name: "SequenceForm",
    data() {
      return {
            read: '',
              db: [],
      mismatches: 5,
           error: '',
         loading: false
      }
    },
    methods: {
      /**
       * performs error checking, then emits form status to App.vue
       */
      submit() { 
        const regex = /[^ATGC]/i;
        if (this.read.search(regex) !== -1) {
          this.error = "only DNA sequences valid";
        } else if (this.read.length < 17) {
          this.error = "input a DNA sequence of length 17 or higher";
        } else if (this.db.length === 0) {
          this.error = "select a database";
        } else {
          this.error = '';
          this.loading = true;
          this.read = this.read.toUpperCase();
          this.$emit("valid");
        }
      }
    },
    emits: ["valid"]
  }
</script>


<template>
  
  <div class="parent">

    <div style="display: flex; justify-content: center; margin: 0;"><a href="http://aspendb.uga.edu/" target="_blank"><img src="../assets/aspendb_bw2.png" style="opacity: 0.9" width="100" height="60  "/></a></div>
    <div class="divider"></div> 
    <p style="width: fit-content; padding: 0 1%; border-radius: 3px; background-color: #f2f2f2">input primer/probe/gRNA sequence</p>
    <textarea rows="2" placeholder="GGGTTCTGCCAATTTAAGCCACATGGCTCAATGGGAGA" v-model="read"></textarea> 
    <div style="display: flex; align-items: center;">
    <p style="width: fit-content; padding: 0 1%; border-radius: 3px; background-color: #f2f2f2">select one or more genomes</p>
    </div>

    <div style="display: flex; flex-direction: column;">
      
            
        <div class="checkboxes">

          <div>
            <input v-model="db" type="checkbox" value="PtrichocarpaV3.1">
            <label>P. trichocarpa V3.1</label>
          </div>
          <div>
            <input v-model="db" type="checkbox" value="PtrichocarpaV4.0">
            <label>P. trichocarpa V4.0</label>
          </div>
          <div>
            <input v-model="db" type="checkbox" value="DeltoidesWV94">
            <label>P. deltoides WV94 V2.1</label>
          </div>
          <div>
            <input v-model="db" type="checkbox" name="717V5" value="717V5">
            <label for="717V5">P. tremula x alba 717 V5</label>
          </div>

      </div>

      <div>
        <p style="width: fit-content; padding: 0 1%; border-radius: 3px; background-color: #f2f2f2">set mismatch number</p>
        <input v-model="mismatches" type="number" class="mismatch">
        <label>max</label>
      </div>

     

      <button v-on:click="submit()">
        <p>search</p>
      </button>

      <!-- display on error -->
      <div class="error" v-if="error != ''">
        <div style="width: fit-content; background-color: #f5f5f7; border-radius: 2px; padding: 0 1%">
          <p> {{ error }} </p>
        </div>
      </div>

      <div style="height: 9px; display: flex; justify-content: center; align-items: center;">
        <div v-if="loading === true" class="lds-ellipsis"><div></div><div></div><div></div><div></div></div>
      </div>

    </div>
    
  </div>
</template>


<style scoped>

p {
    margin: 1% 2% 1% 2%;
    font-size: 10px;
}
label {
    font-size: 10px;
}
textarea {
    width: 96%;
    margin: 0 2% 0 2%;
    padding: 4%;
    font-size: 12px;
    font-weight: 500px;
    border: 1px solid var(--color-background);
    resize: none;
}
input {
    margin: 1% 2%;
    width: fit-content;
}
button {
    width: 45%;
    margin: 4% auto 5% auto;
    color: var(--vt-c-black-mute);
    background-color: #f5f5f7; 
    border-radius: 1px;
    border: 1px solid var(--vt-c-black-mute);
    cursor: pointer;
}
button:hover {
    background-color: #f2f2f2;
}
.parent {     
    width: 90vw;                                                                                                                                             
    display: flex; 
    flex-direction: column;
    flex-shrink: 0;
    margin: 4%;
    background-color: #fbfbfb;
    border: 1px solid var(--color-background);

}
.divider {
    width: 90%; 
    border: 1px solid var(--vt-c-black-mute); 
    border-bottom-style: none; 
    margin: 0 5% 2% 5%;
}
.checkboxes {
    width: 100%;
    margin: 0 2%;
    display: flex;
    flex-direction: column;
}
.mismatch {
  width: 50px;
}
.error {
    margin-top: 1%;
    display: flex;
    justify-content: center;
    align-items: center;
}



.lds-ellipsis {
  display: flex;
  position: relative;
  width: 40px;
  height: 40px;
}
.lds-ellipsis div {
  position: absolute;
  top: 27px;
  width: 7px;
  height: 7px;
  border-radius: 50%;
  background: #222222;
  opacity: 0.85;
  animation-timing-function: cubic-bezier(0, 1, 1, 0);
}
.lds-ellipsis div:nth-child(1) {
  left: -9px;
  animation: lds-ellipsis1 0.7s infinite;
}
.lds-ellipsis div:nth-child(2) {
  left: -9px;
  animation: lds-ellipsis2 0.7s infinite;
}
.lds-ellipsis div:nth-child(3) {
  left: 15px;
  animation: lds-ellipsis2 0.7s infinite;
}
.lds-ellipsis div:nth-child(4) {
  left: 39px;
  animation: lds-ellipsis3 0.7s infinite;
}
@keyframes lds-ellipsis1 {
  0% {
    transform: scale(0);
  }
  100% {
    transform: scale(1);
  }
}
@keyframes lds-ellipsis3 {
  0% {
    transform: scale(1);
  }
  100% {
    transform: scale(0);
  }
}
@keyframes lds-ellipsis2 {
  0% {
    transform: translate(0, 0);
  }
  100% {
    transform: translate(24px, 0);
  }
}



@media (min-width: 1024px) {
.parent {
    width: 70vw; 
    margin: 5% auto 2% auto;
    padding: 1rem 8rem 2rem 8rem;
    opacity: .95;
    border-radius: 3px;
    box-shadow: 0px 0px 1px 0px var(--color-background); 
    border: none;

}
.checkboxes {
    margin: 0;
}
p {
    margin: 1% 0 1% 0;
    font-size: 18px;
}
textarea {
    width: 100%;
    padding: 4px;
    margin: 0;
    background-color: rgba(247, 247, 247, 0.90);
    font-size: 16px;
    font-weight: 500px;
}
input {
    margin: 0.5% 1%;
    width: fit-content;
}
label {
    color: var(--vt-c-black-mute);
    font-size: 16px;
}
button {
    width: 25%;
    padding: 5px;
    margin: 1% 0;
    font-size: 14px;
    margin-left: auto;
    margin-right: auto;
}
.divider {
    width: 100%; 
    margin: 0 0 2% 0; 
} }

@media (min-width: 2048px) {
.parent {
    padding: 4rem 20rem 4rem 20rem;
}
.divider {
    margin-bottom: 0;
}
textarea {
    width: 1000px;
}
button {
    width: 12.5%;
}
input {
    margin: 0.5%;
} }



</style>
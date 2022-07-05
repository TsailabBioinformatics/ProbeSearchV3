<script>
  export default {
    name: "AlignmentResult",
    props: ['data', 'id', 'total'],
    methods: {
      swapleft() {this.$emit("swap-left") },
      swapright() { this.$emit("swap-right") },
      reload() { 
        window.location.reload();
        window.scrollTo(0, 0);
      },
      expand() {
        // emit signal to app, instantiating new component
        // edit current component
        this.$emit("expand")
      },
      download() { // TODO 
      }
    },
    emits: ["swap-left", "swap-right", "expand"]
  }
</script>

<template>
  <div class="parent">

    <p class="expand" @click="expand()"> ⛶ </p>
    <div class="header"> 
        <div>
          <p class="swap-left" v-if="id > 1" @click="swapleft()"> &lt; </p> <p v-else> &nbsp; </p> 
          <p class="previous-instructions"> previous result </p>
        </div>
        <h3 v-if="this.total > 1">Alignment Result {{id}} </h3> 
        <div>
          <p class="swap-right" v-if="id < total" @click="swapright()"> > </p> <p v-else> &nbsp;</p>
          <p class="next-instructions"> next result </p>
        </div>
    </div>

    <pre>
        <p>{{data}}</p>
    </pre>

    <div class="footer">
        <!-- TODO <p class="download" @click="download"> download </p> -->
        <p class="search-again" @click="reload"> search again </p>
    </div>
    
  </div>
</template>

<!-- css -->
<style scoped>
  .parent {
      width: 90vw;  
      display: flex; 
      flex-direction: column;
      flex-shrink: 0;
      margin: 0;
      background-color: var(--vt-c-white-mute);
      /* border: 1px solid hsla(160, 100%, 37%, 1); */
  }
  .header {
    display: flex;
    justify-content: space-evenly;
    align-items: center;
  }
  .swap-left:hover {
    cursor: pointer;
  }
  .swap-right:hover {
    cursor: pointer;
  }
  .footer {
    display: flex;
    align-items: flex-end;
    justify-content: flex-end; /* TODO - add download - space-between */
  }
  .download:hover {
    text-decoration: underline;
    cursor: pointer;
  }
  .search-again {
    padding: 0 1%;
    background-color: #f2f2f2;
    border-radius: 2px;
  }
  .search-again:hover {
    cursor: pointer;
  }
  .previous-instructions {
      display: none;
  }
  .next-instructions {
    display: none;
  }
  .swap-left:hover + .previous-instructions {
      display: flex;
      padding: 5px;
      width: 90px;
      font-size: 8px;
      text-align: center;
      color: rgba(60, 60, 60, 0.66);
      background-color: rgba(247, 247, 247, 0.90);
      border: 1px solid var(--color-background);
      position: absolute;
      left: -50px;
      bottom: 15px;
  }
  .swap-right:hover + .next-instructions {
      display: flex;
      padding: 5px;
      width: 65px;
      font-size: 8px;
      text-align: center;
      color: rgba(60, 60, 60, 0.66);
      background-color: rgba(247, 247, 247, 0.90);
      border: 1px solid var(--color-background);
      position: absolute;
      right: -50px;
      bottom: 15px;
  }
  .expand {
    position: absolute;
    right: 2px;
    top: 2px;
  }
  .expand:hover {
    cursor: pointer;
  }
  p {
      margin: 1% 2%;
      font-family: 'DM Mono', monospace; 
      font-size: 10px;
      font-weight: 200;
  }
  button {
      width: 20%;
      height: 15%;
      color: rgba(60, 60, 60, 0.66);
      background-color: #f2f2f2; 
      border-radius: 5px;
      border: none;
      cursor: pointer;
  }
@media (min-width: 1024px) {
  p {
      margin: 5px;
      font-family: 'DM Mono', monospace; 
      font-size: 16px;
      font-weight: 300;
      color: var(--vt-c-black-mute);
  }
  h3 {
    color: var(--vt-c-black-mute);
  }
  .parent {
      width: 70vw;
      margin: 1% 3% 15% 3%;
      padding: 2rem 4rem 2rem 4rem;
      display: flex; 
      flex-direction: column;
      flex-shrink: 0;
      background-color: #fbfbfb;
      opacity: .95;
      box-shadow: 0px 0px 1px 0px var(--color-background); 
      border: 1px solid var(--vt-c-black-mute);
      border-radius: 2px;
  }
  .swap-left:hover + .previous-instructions {
      display: flex;
      padding: 5px;
      width: 105px;
      font-size: 10px;
      text-align: center;
      color: rgba(60, 60, 60, 0.66);
      background-color: rgba(247, 247, 247, 0.90);
      border: 1px solid var(--color-background);
      position: absolute;
      left: -75px;
      bottom: 25px;
  }
  .swap-right:hover + .next-instructions {
      display: flex;
      padding: 5px;
      width: 80px;
      font-size: 10px;
      text-align: center;
      color: rgba(60, 60, 60, 0.66);
      background-color: rgba(247, 247, 247, 0.90);
      border: 1px solid var(--color-background);
      position: absolute;
      right: -50px;
      bottom: 25px;
  }

}
</style>
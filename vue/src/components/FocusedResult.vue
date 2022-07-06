<script>
  export default {
    name: "FocusedResult",
    props: ['data', 'id'],
    methods: {
      contract() {
        this.$emit("contract")
      },
      download() { /* TODO */ }
    },
    emits: ["contract"]
  }
</script>


<template>

  <div class="parent">

        <p class="contract" @click="contract()"> ⛶ </p>
        <p class="instruction"> minimize </p>
        <pre>
            <p>{{data}}</p>
        </pre>
    
  </div>

</template>


<style scoped>

.parent {
    width: 90vw;
    height: 100%;
    min-height: 90vh;
    padding: 2rem 4rem 2rem 4rem;
    background-color: #fbfbfb;
    box-shadow: 0px 0px 1px 0px var(--color-background); 
    border: 1px solid var(--vt-c-black-mute);
    border-radius: 2px;
}
.contract {
    position: absolute;
    right: 6px;
    top: 6px;
}
.contract:hover {
    cursor: pointer;
}
.instruction {
    display: none;
}

@media (min-width: 1024px) {
.parent {
    padding: 2rem 4rem 2rem 4rem;
    margin: 5vh 5vw 5vh 5vw;
} 
.contract:hover + .instruction {
    display: flex;
    padding: 5px;
    width: 60px;
    font-size: 10px;
    text-align: center;
    color: rgba(60, 60, 60, 0.66);
    background-color: rgba(247, 247, 247, 0.90);
    border: 1px solid var(--color-background);
    position: absolute;
    right: -40px;
    top: -25px;
} }

</style>
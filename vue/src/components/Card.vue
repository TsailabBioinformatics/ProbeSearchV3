<script>
  export default {
    name: "Card",
    props: ['title', 'id', 'header', 'result'],
    methods: {
      maximize() { this.$emit("fullscreen") }
    },
    emits: ["fullscreen"]
  }
</script>


<template>

  

  <div class="card" @click="maximize()">

        <h1> {{ title }} </h1>
    
  </div>

</template>


<style scoped>

.card {
    display: flex;
    justify-content: center;
    align-items: center;
    width: 225px;
    height: 275px;
    margin: 0 1% 1% 1%;
    padding: 2%;
    background-color: #fbfbfb;
    opacity: .95;
    box-shadow: 0px 0px 1px 0px var(--color-background); 
    border-radius: 3px;
}
.card:hover {
  background-color: #f2f2f2;
  cursor: pointer;
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
p, h3 {
  font-family: 'Oswald', sans-serif;
}

@media (min-width: 1024px) {
.parent {
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
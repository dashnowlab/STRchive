import React, { useEffect, useRef } from 'react';
import * as echarts from 'echarts';
import data from "../../../data/plots/gnomad-data.json"

const GnomADChart = ({id, sex}) => {
  const chartRef = useRef(null);

  useEffect(() => {
    if (chartRef.current) {
      const myChart = echarts.init(chartRef.current);
      var plotdata;
      if (sex == undefined) { plotdata = data[id]; }
      else { plotdata = data[id][sex]; }

      const means = plotdata.xvalues;
      const lowerbounds = plotdata.confidence_lowerbounds;
      const upperbounds = plotdata.confidence_upperbounds;
      const ylabels = plotdata.ylabels;
      
      const confvalues = lowerbounds.map((x, i) => {
        return upperbounds[i] - lowerbounds[i];
      });

      const option = {
        grid: {
          containLabel: true,
        },
        title: {
          text: plotdata.title,
        },
        yAxis: {
          type: 'category',
          data: ylabels,
          margin: 100,
          axisLabel: {
            fontSize: 15,
          }  
        },
        toolbox: {
          // feature: {
          //   saveAsImage: { type: 'png' },
          //   restore: {},
          //   dataZoom: {},
          // }
        },
        tooltip: {
          trigger: 'axis',
          axisPointer: {
            type: 'shadow'
          },
          formatter: function (params) {
            var tar = params[0];
            var text = tar.name + '<br/>' + tar.seriesName + ' : ' + tar.value.toFixed(3);
            text += '<br/>95% Confidence Interval';
            text += '<br/>Lower Bound: ' + params[1].value.toFixed(3);
            text += '<br/>Upper Bound: ' + (params[1].value+params[2].value).toFixed(3);
            return text;
          }
        },
        xAxis: {
          type: 'value',
          name: 'Pathogenic Genotype (%)',
          nameLocation: 'centre',
          nameGap: 15,
          nameTextStyle: {
            fontSize: 16,
            verticalAlign: 'top',
            padding: [30, 0, 0, -65]
          },
          offset: 0,
        },
        series: [
          {
            name: "Percentage Pathogenic",
            data: means,
            type: 'scatter'
          },
          {
            data: lowerbounds,
            barWidth: 0,
            type: 'bar',
            stack: 'Total',
            itemStyle: {
              borderColor: 'transparent',
              color: 'transparent'
            },
            emphasis: {
              itemStyle: {
                borderColor: 'transparent',
                color: 'transparent'
              }
            },
          },
          {
            data: confvalues,
            barWidth: 3,
            type: 'bar',
            stack: 'Total',
          },
        ]
      };

      myChart.setOption(option);
    }
  }, []);

  return <div ref={chartRef} style={{ width: '600px', height: '400px' }} />;
};

export default GnomADChart;
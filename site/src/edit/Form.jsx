import JsonForm from "@rjsf/core";
import validator from "@rjsf/validator-ajv8";
import schema from "~/STRchive-loci.schema.json";

const Form = ({ data }) => {
  return (
    <JsonForm
      schema={schema}
      validator={validator}
      formData={data}
      onChange={console.info}
    />
  );
};

export default Form;

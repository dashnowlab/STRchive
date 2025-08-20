import { Fragment, useCallback, useState } from "react";
import { LuSend } from "react-icons/lu";
import { groupBy, mapValues } from "lodash-es";
import Button from "@/components/Button";
import FormWrapper from "@/components/Form";
import Heading from "@/components/Heading";
import Help from "@/components/Help";
import TextBox from "@/components/TextBox";
import { validate } from "@/data/schema";
import classes from "./Form.module.css";
import {
  properties,
  required as requiredKeys,
} from "~/STRchive-loci.schema.json";

const sections = groupBy(
  Object.entries(properties).map(([key, value]) => ({ key, ...value })),
  "section",
);

const Form = ({ data: initialData = {} }) => {
  const [data, _setData] = useState(initialData);
  const [errors, _setErrors] = useState(null);

  const setErrors = useCallback(
    (newData) => {
      let errors = validate(newData ?? data) ?? [];
      errors = errors.map(({ instancePath, message, params }) => [
        instancePath.split("/").pop() || params.missingProperty,
        message,
      ]);
      errors = Object.fromEntries(errors);
      _setErrors(errors);
    },
    [data],
  );

  const setData = useCallback((key, value) => {
    _setData((data) => {
      const newData = { ...data, [key]: value };
      setErrors(newData);
      return newData;
    });
  }, []);

  return (
    <FormWrapper onSubmit={console.info}>
      <div></div>
      {Object.entries(sections).map(([section, properties]) => {
        if (section === "undefined") return <Fragment key={section} />;
        return (
          <section key={section}>
            <Heading level={2}>{section}</Heading>

            <div className={classes.form}>
              {properties.map((property) => {
                const {
                  key,
                  type,
                  title,
                  description,
                  enum: _enum,
                  examples,
                  pattern,
                } = property;

                const types = [type].flat();

                const allowNull = types.includes("null");

                const required = requiredKeys.includes(key) && !allowNull;

                const tooltip = [
                  description,
                  ...(examples?.length ? [" "] : []),
                  examples?.length ? `Examples: ${examples.join(", ")}` : "",
                ]
                  .filter(Boolean)
                  .join("<br/>");
                const label = (
                  <>
                    {title ?? key}
                    {description && <Help>{tooltip}</Help>}
                  </>
                );

                const errorMessage = errors?.[key] || "";
                const error = !!errorMessage;
                const errorRow = error ? (
                  <span className={classes.error}>{errorMessage}</span>
                ) : null;

                let field = <></>;

                if (!_enum && types.includes("string"))
                  field = (
                    <TextBox
                      key={key}
                      ref={(el) => {
                        if (error) el?.setCustomValidity(errorMessage);
                        else el?.setCustomValidity("");
                      }}
                      name={key}
                      label={label}
                      value={data[key] ?? ""}
                      pattern={pattern}
                      required={required}
                      error={errors?.[key]}
                      onChange={(value) => {
                        if (value === "" && allowNull) setData(key, null);
                        else if (value === "" && !allowNull)
                          setData(key, undefined);
                        else setData(key, value);
                      }}
                    />
                  );
                else return <Fragment key={key} />;

                return (
                  <Fragment key={key}>
                    {field}
                    {errorRow}
                  </Fragment>
                );
              })}
            </div>
          </section>
        );
      })}

      <section>
        <Button type="submit" design="bubble" onClick={() => setErrors()}>
          <LuSend />
          <span>Submit</span>
        </Button>
      </section>
    </FormWrapper>
  );
};

export default Form;
